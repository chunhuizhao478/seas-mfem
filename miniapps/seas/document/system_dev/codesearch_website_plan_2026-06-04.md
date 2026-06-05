# Local Code-Search + Q&A Website for the SEAS / MFEM Codebase — Plan

- **Date:** 2026-06-04
- **Author:** Chunhui Zhao (plan drafted with Claude Code)
- **Status:** Plan / not yet implemented
- **Scope decision:** Hybrid (lexical code-nav + semantic NL Q&A), single-user, runs locally on the Mac.
- **Lives outside the MFEM repo** (a sibling directory) so it never touches the build or the git history.

---

## 1. Goal

Stand up a personal website that lets me both:

1. **Navigate the code** — fast regex / substring / symbol search and jump-to-definition over
   the full C++/Python tree (lexical, deterministic).
2. **Ask questions in English** — "why must the friction solver use Brent not Newton?",
   "where is the fault-locality partition built?" — and get a synthesized answer with
   citations into the source/docs (semantic RAG).

The two are different engines under the hood; this plan unifies them behind one small
web app so there is a single search box and a single "ask" box.

---

## 2. Codebase profile (drives every sizing decision)

Measured on branch `system/spatial_dyn_driver`, 2026-06-04 (`git ls-files`):

| Category | Count / size | Note |
|---|---|---|
| C++ source (`.cpp .hpp .h .inl .cc`) | ~1,460 files, ~878K LOC (with Python) | **Mostly the full MFEM library**; my code is `miniapps/seas/**` |
| Markdown docs (`.md`) | 410 files, ~209K LOC | Debug history, porting reports, BP5-debug v1–v62, runbooks — **the crown jewel for semantic search** |
| Python (`.py`) | 85 files | Post-processing, mesh, plotting, estimators |
| Job scripts (`.sbatch`) | 337 | Frontera run recipes |
| Config (`.toml`) | 37 | Benchmark/scheme configs |
| Data artifacts (`.vtu .pvtu .vtk .png .pdf .mesh .msh .dat .csv .geo .step`) | ~1,400 files | **Must be excluded** — ~90% of the 2.3 GB tree, pure noise in code results |
| Tracked files total | 4,103 | |
| Working tree | 2.3 GB | |
| `.git` | 1.2 GB | |

**Key implications:**

- The index must be **scoped and filtered**, or it bloats with binary artifacts and the
  results are polluted.
- The repo is the full MFEM library **plus** the seas miniapp. Index both (MFEM is useful
  reference) but **tag and rank `miniapps/seas/**` first**.
- The 209K-LOC markdown corpus is exactly the content where NL Q&A beats grep — it encodes
  tribal knowledge (sign conventions, why-Brent-not-Newton, the σ_n drift history) that a
  literal search cannot surface.

---

## 3. Architecture decision

**Chosen:** one lightweight **FastAPI** app that serves both a lexical search box and a
semantic "ask" box. No heavyweight code-search server (Sourcegraph / OpenGrok) because this
is single-user and local — those shine for teams and add operational weight we do not need.

### Component stack

| Concern | Tool | Why (local / Apple Silicon / single-user) |
|---|---|---|
| Lexical search | **ripgrep** (`rg --json`) as a subprocess | No index to maintain; searches the 1M-LOC tree in < 1 s |
| Symbol jump-to-def | **universal-ctags** → `tags` file | Cheap C++/Python definitions without a full clangd/LSP setup |
| Chunking for embeddings | **tree-sitter** (`tree-sitter-languages`) | C++/Py split by function; markdown split by heading |
| Embeddings | **Ollama** `nomic-embed-text` | Native arm64, free, **fully local** (nothing leaves the machine at index time) |
| Vector store | **LanceDB** | Embedded (a file on disk), no server process to run |
| Hybrid ranking | reciprocal-rank fusion of BM25 (ripgrep hit) + vector similarity | ~20 lines, no extra infra |
| Answer / synthesis | **Claude API** (`anthropic`) | Cited NL answers; lower friction than wiring a local LLM |
| Backend | FastAPI + uvicorn | Matches the project's Python stack |
| Frontend | single `index.html`, syntax highlight via **Shiki** (or Monaco) | Minimal; line-anchor permalinks |
| Deploy | `uvicorn app:app` on `localhost:8000` | Personal; no auth, no container required |

### Privacy model (matters because this is "just me")

- With **local Ollama embeddings**, *indexing is fully local* — no code is sent anywhere.
- Only the final `/ask` call sends the **retrieved snippets** (top-k chunks) to the Claude API
  for synthesis.
- If even that must stay offline, swap the answer step for a local Ollama LLM
  (e.g. `qwen2.5-coder`) — lower answer quality but zero data egress.

---

## 4. Directory layout (sibling dir, NOT inside the MFEM repo)

```
~/projects/seas-codesearch/
  config.toml      # repo path, include/exclude globs, seas-priority boost, model names
  index.py         # git ls-files -> filter -> chunk -> embed -> LanceDB; also runs ctags
  app.py           # FastAPI: /search /symbol /ask /file
  web/index.html   # two boxes (search + ask), results pane, citation deep-links
  .venv/
  data/            # LanceDB table + tags file (gitignored / local only)
```

Keeping it a sibling directory means it never interferes with the MFEM build
(`compile_commands.json`, `make`) and never gets committed into the research repo.

---

## 5. Bootstrap

```bash
brew install ripgrep universal-ctags ollama
ollama pull nomic-embed-text
cd ~/projects/seas-codesearch
python3 -m venv .venv && source .venv/bin/activate
pip install fastapi uvicorn lancedb anthropic tree-sitter-languages tomli
export ANTHROPIC_API_KEY=...        # only needed for the /ask synthesis step
```

---

## 6. Indexing scope and filters (the part specific to THIS repo)

This is the make-or-break configuration. Put it in `config.toml`.

- **Repo root:** `/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver`
- **Include extensions:** `cpp hpp h inl cc py md toml sbatch cmake`
- **Exclude extensions:** `vtu pvtu vtk png pdf mesh msh dat csv geo step stp xyz`
- **Exclude path globs:** `**/gold/**`, `**/out_*/**`, `**/results/**`, `**/plots/**`,
  `**/*.prephase*bak`, `**/.git/**`
- **Tag** each chunk `seas` if the path starts with `miniapps/seas/`, else `mfem`;
  **boost `seas` chunks** in both lexical and vector ranking.
- **RAG priority corpus** (index first, weight up): `miniapps/seas/document/**`,
  `miniapps/seas/debug_document/**`, and all `CLAUDE.md` / `ARCHITECTURE.md` files — this is
  the highest-value text for question answering.

Enumerate the file list with `git ls-files` (respects tracking, skips build artifacts and
the gitignored `.msh`/`out_*`), then apply the include/exclude filter on top.

---

## 7. Components in detail

### 7.1 `index.py`
1. `git ls-files` at the repo root → apply include/exclude filter.
2. For each file:
   - **code** (`.cpp/.hpp/.h/.inl/.cc/.py`): tree-sitter → one chunk per function/method/class,
     carrying `{path, start_line, end_line, lang, scope=seas|mfem, symbol_name}`.
   - **markdown** (`.md`): split by heading (H1–H3) → one chunk per section, carrying
     `{path, heading_path, start_line, end_line, scope}`.
3. Embed each chunk via Ollama `nomic-embed-text` → upsert into a LanceDB table with metadata.
4. Run `ctags -R --languages=C,C++,Python -f data/tags <included files>` for the symbol map.
5. Write an index manifest `{indexed_at, n_chunks, n_files, git_commit}`.

Rough cost on this tree: a few minutes (embeddings dominate; ~tens of thousands of chunks).

### 7.2 `app.py` (FastAPI endpoints)
- `GET /search?q=&regex=&scope=` → shell out to `rg --json` (scoped to include globs,
  excludes applied) → return `[{path, line, preview}]`. Instant; no index touched.
- `GET /symbol?name=` → look up in `data/tags` → return definition site(s) for jump-to-def.
- `GET /ask?q=` → **hybrid retrieve** (see §8) → top-k chunks → Claude API with a
  "answer ONLY from these snippets, cite path:line" system prompt → stream answer + citations.
- `GET /file?path=&line=` → read file, syntax-highlight with Shiki, return HTML anchored at
  the requested line (the target of every citation and search-result click).

### 7.3 `web/index.html`
- Top: a **Search** box (lexical) and an **Ask** box (RAG), toggle or side-by-side.
- Search results: a table (path:line + preview), click → `/file` viewer at the line.
- Ask answer: streamed text + **citation chips** that deep-link into the `/file` viewer.
- A `scope` toggle: "seas only" vs "all (incl. MFEM)".

---

## 8. Hybrid retrieval (for `/ask`)

For a question `q`:
1. **Lexical candidate set:** run `rg` for the salient terms → collect the hit files/lines,
   map to overlapping chunks (BM25-ish lexical score).
2. **Vector candidate set:** embed `q` via Ollama → LanceDB top-N by cosine.
3. **Fuse** with reciprocal-rank fusion: `score(c) = Σ_lists 1 / (k + rank_list(c))`, `k≈60`,
   with a multiplicative boost for `scope == seas`.
4. Take top-k fused chunks (e.g. 8–12), pack into the Claude prompt with their `path:line`
   labels, require citations.

This keeps exact-term questions (symbol names, flags) sharp via lexical, and conceptual
questions ("why does σ_n drift") strong via vectors.

---

## 9. Build phases

1. **Lexical half first** (fastest to something usable): `/search` + `/symbol` + the `/file`
   viewer + the UI. No embeddings needed; immediately useful for day-to-day navigation.
2. **Indexing pipeline:** `index.py` (chunk + embed + LanceDB + ctags).
3. **Semantic half:** `/ask` with hybrid retrieval + Claude synthesis; wire citation chips.
4. **Refresh strategy:** re-run `index.py` by hand, or add a `git post-commit` hook in the
   MFEM repo that re-indexes incrementally (only changed files since last manifest commit).

---

## 10. Refresh / maintenance

- Lexical (`/search`, `/symbol`) needs no rebuild for `rg`; regenerate `tags` after large edits.
- Semantic needs re-embedding of changed chunks. Track `git_commit` in the manifest and only
  re-embed files whose blob hash changed since the last index.
- Everything is local files (LanceDB table + tags); nuke `data/` to rebuild from scratch.

---

## 11. Future extensions

- **Real C++ cross-references** ("find all references", not just definitions): generate
  `compile_commands.json` from the MFEM build and run **`scip-clang`** → load SCIP for
  precise xref. This is the hardest piece; deferred until plain ctags proves insufficient.
- **Sharing with collaborators:** move from localhost to a small VM behind Tailscale or
  basic-auth; add a re-index-on-push GitHub Action. (Out of scope for the single-user goal.)
- **Local-only answers:** replace the Claude synthesis step with a local Ollama LLM for full
  offline operation.
- **Knowledge-graph view:** the `understand-anything:understand` agent already builds a
  browsable graph of the codebase — a complementary navigation mode if graph browsing is
  wanted over a search box.

---

## 12. Off-the-shelf alternatives considered (and why not, for now)

| Option | What it gives | Why deferred |
|---|---|---|
| **Sourcebot** (self-host, Zoekt-based, +AI) | Closest single product to "a code-search website" | Geared to indexing remote git hosts; more than needed for one local repo |
| **Sourcegraph** (self-host) | Search + jump-to-def + blame, polished | Heaviest; OSS edition recently relicensed/trimmed; team-oriented |
| **OpenGrok** | Free, mature C++ symbol cross-refs | Tomcat/Java footprint; no semantic Q&A |
| **Hound / livegrep** | Instant regex web UI in a single Go binary | Lexical only; afternoon MVP but no RAG, no unified UI |
| **Bloop** | Semantic code search with local embeddings | Less maintained now |

The custom FastAPI app is preferred here because it is the only path that puts **both**
lexical and semantic search behind **one** local UI in ~250–350 lines, with no server to
operate and full control over the seas-vs-MFEM scoping.

---

## 13. Open decisions

- Index all of MFEM, or `miniapps/seas/**` only? (Plan: index both, rank seas first.)
- Answer model: Claude API (default) vs local Ollama LLM (offline)?
- Chunk granularity for C++ headers full of templates (function-level may be too fine for
  MFEM core) — may need a fallback to file/section chunking for very dense headers.

---

## 14. One-line summary

A ~300-line local FastAPI app: **ripgrep + ctags** for instant lexical code-nav and
**tree-sitter → Ollama embeddings → LanceDB → Claude** for cited NL Q&A, behind one web UI,
with the data-artifact filter and `miniapps/seas/**` ranking boost that this repo specifically
needs.
