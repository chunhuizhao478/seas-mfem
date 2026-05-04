# CGAL header shim for cgal-cpp 5.6.1

## Why this directory exists

`conda-forge::cgal-cpp 5.6.1` (the package the laptop recipe in
`tools/README.md` installs) ships a header
`<CGAL/boost/graph/iterator.h>` whose `operator bool()` overloads
on three classes —

- `Halfedge_around_source_iterator` (line 217),
- `Halfedge_around_target_iterator` (line 311),
- `Halfedge_around_face_iterator`   (line 403)

— call `this->base()`. Those classes do **not** inherit from
`boost::iterator_facade`; they have only the data members
`{anchor, pos, g, winding}`. There is therefore no `base()`
member and Apple Clang fails the build with:

```
error: no member named 'base' in
        'Halfedge_around_target_iterator<Graph>'
```

This error fires the moment any code in the corefinement /
PMP include chain instantiates one of these iterators.

## What the shim does

A copy of `iterator.h` lives at
`tools/cgal_patch/CGAL/boost/graph/iterator.h`. The three offending
bodies are rewritten to use the iterator's actual sentinel field:

```diff
-     return (! (this->base() == nullptr));
+     return (this->g != nullptr);
```

Default-constructed iterators set `g(nullptr)` (see line 287 etc. of
the original header), so `this->g != nullptr` reproduces the intent
of the original `operator bool()`.

`tools/CMakeLists.txt` puts this shim at the FRONT of the include
path with `target_include_directories(... BEFORE ...)`, so the patched
header shadows the broken one for the `corefine_faults_hello` (and
future Phase 1 `corefine_faults`) targets only. Other binaries
that link against the conda env's CGAL are unaffected.

## When to remove this shim

Delete `tools/cgal_patch/` and the
`target_include_directories(... BEFORE ...)` line in
`tools/CMakeLists.txt` once **either** of the following holds:

- conda-forge ships a fixed `cgal-cpp` 5.6.x or 5.7.x package
  (verify by `grep "this->base() == nullptr"` returning empty in
  `$CONDA_PREFIX/include/CGAL/boost/graph/iterator.h`); or
- the project moves to CGAL ≥ 6.0 (out of the currently-tested
  ABI window — would also require updating
  `find_package(CGAL 5.4)` and the version guard in
  `tools/CMakeLists.txt`).

Until then this shim is the simplest way to keep the conda env
clean while still compiling the CGAL corefinement headers.

## Provenance

The original header was copied verbatim from
`$CONDA_PREFIX/include/CGAL/boost/graph/iterator.h` of the
`cgal-cpp 5.6.1 h6a79a76_1` conda-forge package on 2026-04-29.
The diff against the original is exactly the three sed
substitutions documented above. No other lines of the file have
been changed; if conda updates the package, re-copy the file and
re-apply the same sed.
