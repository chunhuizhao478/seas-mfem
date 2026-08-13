#!/usr/bin/env python3
"""medit_hdr.py -- locate the sections of a MEDIT .mesh file instead of hardcoding them.

Why this exists: the tools vendored from the PREFERRED tree carried that mesh's
byte layout as constants --

    LV, LT, LS = 7, 13428231, 84344339
    NV, NT, NS = 12042113, 70641218, 3049525

On any other mesh those line numbers land in the middle of the vertex block, and
because `pandas.read_csv(skiprows=..., nrows=...)` will happily parse whatever is
there, the failure is SILENT: coordinates read as connectivity, no exception, a
plausible-looking mesh out the far end. That is the worst class of bug in this
pipeline, so the offsets are measured from the file itself.

Fast on multi-GB files: mmap, find the keyword's byte offset, then count newlines
up to it with numpy rather than iterating lines in Python.

Returns {name: (line_of_count, count)} where the first DATA row is line+1, matching
the `skiprows=line+1, nrows=count` convention the callers already use.
"""
import mmap
import re

import numpy as np

SECTIONS = ("Vertices", "Triangles", "Tetrahedra")


def medit_sections(path, names=SECTIONS):
    out = {}
    with open(path, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        buf = np.frombuffer(mm, dtype=np.uint8)
        for name in names:
            key = ("\n" + name + "\n").encode()
            pos = mm.find(key)
            if pos < 0:                       # some writers pad with \r or spaces
                m = re.search((r"\n\s*" + name + r"\s*\n").encode(), mm)
                if m is None:
                    continue
                pos = m.start()
            # line index of the keyword = number of newlines strictly before it,
            # counting the one at `pos` as terminating the PREVIOUS line
            kw_line = int(np.count_nonzero(buf[:pos + 1] == 0x0A))
            # kw_line is the 0-based line of the KEYWORD. Callers use
            # `skiprows=L+1, nrows=n`, so L must be the COUNT line, one below it.
            # Validated against PREFERRED's hardcoded LV/LT/LS = 7 / 13,428,231 /
            # 84,344,339: without this +1 every section reads one row early.
            cnt_line = kw_line + 1
            nl = mm.find(b"\n", pos + len(key))
            count = int(mm[pos + len(key):nl].split()[0])
            out[name] = (cnt_line, count)
        # `buf` is a zero-copy view onto the mmap, so it must go before close()
        # or mmap raises BufferError("cannot close exported pointers exist").
        del buf
        mm.close()
    return out


if __name__ == "__main__":
    import sys
    for k, (ln, n) in medit_sections(sys.argv[1]).items():
        print(f"{k:12s} count line {ln:>12,}   n = {n:,}   first data row {ln+1:,}")
