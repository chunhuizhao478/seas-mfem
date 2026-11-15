// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// toml_reader.hpp — strict TOML helpers shared by the pff case parsers.
//
// toml11 v3 (extern/toml11).  Reals accept TOML integers too (toml11's
// find_or<double> would otherwise fall back to the default on `E = 50`).
// Check() throws on keys a table does not allow.

#ifndef MFEM_PFF_TOML_READER_HPP
#define MFEM_PFF_TOML_READER_HPP

#include "mfem.hpp"
#include "toml.hpp"

#include <array>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace mfem
{
namespace pff
{

struct TomlReader
{
   using Keys = std::set<std::string>;

   static void Check(const toml::value &t, const std::string &table, const Keys &allowed)
   {
      for (const auto &kv : toml::get<toml::table>(t))
      {
         if (!allowed.count(kv.first))
         {
            throw std::invalid_argument("unknown key '" + kv.first + "' in [" + table + "]");
         }
      }
   }

   static real_t Real(const toml::value &t, const std::string &k, real_t def)
   {
      if (!t.contains(k)) { return def; }
      const toml::value &v = t.at(k);
      if (v.is_integer()) { return static_cast<real_t>(v.as_integer()); }
      return static_cast<real_t>(toml::get<double>(v));
   }

   static int Int(const toml::value &t, const std::string &k, int def)
   {
      return t.contains(k) ? static_cast<int>(toml::get<toml::integer>(t.at(k))) : def;
   }

   static bool Bool(const toml::value &t, const std::string &k, bool def)
   {
      return t.contains(k) ? toml::get<bool>(t.at(k)) : def;
   }

   static std::string Str(const toml::value &t, const std::string &k, const std::string &def)
   {
      return t.contains(k) ? toml::get<std::string>(t.at(k)) : def;
   }

   /// A list of boundary names or attribute numbers, kept as strings.
   static std::vector<std::string> Names(const toml::value &t, const std::string &k)
   {
      std::vector<std::string> out;
      if (!t.contains(k)) { return out; }
      for (const auto &v : toml::get<toml::array>(t.at(k)))
      {
         out.push_back(v.is_integer() ? std::to_string(v.as_integer())
                       : toml::get<std::string>(v));
      }
      return out;
   }

   static std::array<real_t, 3> Vec3(const toml::value &t, const std::string &k,
                                     std::array<real_t, 3> def)
   {
      if (!t.contains(k)) { return def; }
      const auto &a = toml::get<toml::array>(t.at(k));
      if (a.size() < 1 || a.size() > 3)
      {
         throw std::invalid_argument("'" + k + "' needs 1 to 3 components");
      }
      std::array<real_t, 3> out{{0.0, 0.0, 0.0}};
      for (size_t i = 0; i < a.size(); i++)
      {
         out[i] = a[i].is_integer() ? static_cast<real_t>(a[i].as_integer())
                                    : static_cast<real_t>(toml::get<double>(a[i]));
      }
      return out;
   }

};

} // namespace pff
} // namespace mfem

#endif
