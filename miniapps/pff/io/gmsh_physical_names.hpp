// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// gmsh_physical_names.hpp — map Gmsh physical-group names to MFEM attributes.
//
// MFEM's Gmsh reader sets element/boundary attributes from the physical tags
// and ignores $PhysicalNames.  MOOSE addresses sidesets by name (spe_hftc:
// borehole, xmin, ..., zmax) or by number (cmame_revision: boundary = 3).  The
// case files keep MOOSE's choice; this reads the $PhysicalNames block (ASCII in
// both ASCII and binary .msh files, versions 2.2 and 4.1) so a name resolves to
// the tag MFEM uses as the attribute.

#ifndef MFEM_PFF_GMSH_PHYSICAL_NAMES_HPP
#define MFEM_PFF_GMSH_PHYSICAL_NAMES_HPP

#include <cctype>
#include <fstream>
#include <istream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mfem
{
namespace pff
{

/// (dimension, name) -> physical tag.
using PhysicalNameMap = std::map<std::pair<int, std::string>, int>;

inline PhysicalNameMap ReadPhysicalNames(std::istream &in)
{
   PhysicalNameMap m;
   std::string line;
   while (std::getline(in, line))
   {
      if (line.rfind("$PhysicalNames", 0) != 0) { continue; }
      int n = 0;
      in >> n;
      for (int i = 0; i < n; i++)
      {
         int dim = 0, tag = 0;
         in >> dim >> tag;
         std::getline(in, line);
         const auto q0 = line.find('"'), q1 = line.rfind('"');
         if (q0 == std::string::npos || q1 <= q0)
         {
            throw std::runtime_error("malformed $PhysicalNames entry");
         }
         m[{dim, line.substr(q0 + 1, q1 - q0 - 1)}] = tag;
      }
      break;
   }
   return m;
}

inline PhysicalNameMap ReadPhysicalNames(const std::string &msh_file)
{
   std::ifstream in(msh_file, std::ios::binary);
   if (!in) { throw std::runtime_error("cannot open mesh file " + msh_file); }
   return ReadPhysicalNames(in);
}

/// Resolve boundary entries (names or numbers as strings) of dimension
/// bdr_dim = dim - 1 to attribute numbers.
inline std::vector<int> ResolveAttributes(const std::vector<std::string> &entries,
                                          const PhysicalNameMap &names, int bdr_dim)
{
   std::vector<int> out;
   for (const std::string &e : entries)
   {
      bool numeric = !e.empty();
      for (char ch : e) { numeric = numeric && std::isdigit(static_cast<unsigned char>(ch)); }
      if (numeric) { out.push_back(std::stoi(e)); continue; }
      auto it = names.find({bdr_dim, e});
      if (it == names.end())
      {
         throw std::invalid_argument("boundary '" + e + "' is not a physical group of dimension " +
                                     std::to_string(bdr_dim) + " in the mesh file");
      }
      out.push_back(it->second);
   }
   return out;
}

} // namespace pff
} // namespace mfem

#endif
