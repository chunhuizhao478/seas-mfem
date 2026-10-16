// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// build_info.hpp — fully dynamic BP5 cycle (plan rev 3, rule M10): the commit
// the binary was built from.
//
// cycle/cycle.mk passes SEAS_CYCLE_COMMIT = the full `git rev-parse HEAD`,
// with "-dirty" appended when a tracked file differs from that commit.  A
// build that does not go through cycle.mk reports "unknown".  The driver
// prints it in its log; from P1 on it also goes into every output header and
// every checkpoint, and a restart refuses a checkpoint written by another
// commit.

#ifndef MFEM_SEAS_CYCLE_BUILD_INFO_HPP
#define MFEM_SEAS_CYCLE_BUILD_INFO_HPP

#include <cstring>

#ifndef SEAS_CYCLE_COMMIT
#define SEAS_CYCLE_COMMIT "unknown"
#endif

namespace mfem
{
namespace seas
{
namespace cycle
{

/// The commit this binary was built from ("<sha>", "<sha>-dirty" or "unknown").
inline const char *BuildCommit() { return SEAS_CYCLE_COMMIT; }

/// True when the binary was built through cycle.mk from a clean tracked tree,
/// so BuildCommit() names exactly the sources that were compiled.
inline bool BuildIsClean()
{
   const char *c = BuildCommit();
   const std::size_t n = std::strlen(c);
   const char *dirty = "-dirty";
   const std::size_t nd = std::strlen(dirty);
   if (std::strcmp(c, "unknown") == 0) { return false; }
   return !(n >= nd && std::strcmp(c + n - nd, dirty) == 0);
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_BUILD_INFO_HPP
