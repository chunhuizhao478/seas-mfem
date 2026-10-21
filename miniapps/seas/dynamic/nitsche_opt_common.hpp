// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_opt_common.hpp — deck helpers shared by the standalone options of
// the Nitsche driver for complex fault geometries (TPV12, TPV15, TPV22;
// complex-geometries plan rev 2, Section 5.1, Chunhui Zhao 2026-10-20).
//
// Every option is a top-level deck table of its own, read only by the
// Nitsche driver (the SAFS parser never sees it).  A table refuses every key
// it does not know, names the key and lists the known ones: the shared SAFS
// sub-parsers skip unknown keys, so an option key misspelt inside them
// would be ignored silently.  The names carry the Opt prefix so they never
// collide with the driver's own helpers.

#ifndef MFEM_SEAS_NITSCHE_OPT_COMMON_HPP
#define MFEM_SEAS_NITSCHE_OPT_COMMON_HPP

#include "mfem.hpp"

#include <toml.hpp>

#include <initializer_list>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace nitsche_opt
{

/// Abort unless `t` is a table whose keys are all in `known`.
inline void OptRefuseUnknownKeys(const toml::value &t, const char *table,
                                 std::initializer_list<const char *> known)
{
   MFEM_VERIFY(t.is_table(), "[" << table << "] must be a table");
   for (const auto &kv : t.as_table())
   {
      bool ok = false;
      for (const char *k : known) { ok = ok || kv.first == k; }
      if (!ok)
      {
         std::ostringstream ks;
         for (const char *k : known) { ks << " " << k; }
         MFEM_ABORT("[" << table << "]: unknown key '" << kv.first << "' (known:"
                    << ks.str() << ")");
      }
   }
}

/// A number (integer or float) at key `key` of table `t`; required.
inline real_t OptReal(const toml::value &t, const char *table, const char *key)
{
   MFEM_VERIFY(t.contains(key), "[" << table << "] needs " << key << " (no default)");
   const toml::value &v = t.at(key);
   if (v.is_floating()) { return static_cast<real_t>(v.as_floating()); }
   if (v.is_integer()) { return static_cast<real_t>(v.as_integer()); }
   MFEM_ABORT("[" << table << "] " << key << " must be a number");
   return 0.0;
}

/// A string at key `key` of table `t`; required.
inline std::string OptStr(const toml::value &t, const char *table, const char *key)
{
   MFEM_VERIFY(t.contains(key), "[" << table << "] needs " << key << " (no default)");
   MFEM_VERIFY(t.at(key).is_string(), "[" << table << "] " << key << " must be a string");
   return t.at(key).as_string();
}

/// A boolean at key `key` of table `t`; required.
inline bool OptBool(const toml::value &t, const char *table, const char *key)
{
   MFEM_VERIFY(t.contains(key), "[" << table << "] needs " << key << " (no default)");
   MFEM_VERIFY(t.at(key).is_boolean(), "[" << table << "] " << key << " must be true or false");
   return t.at(key).as_boolean();
}

/// An array of numbers at key `key` of table `t`; required, may be empty.
inline std::vector<real_t> OptRealArray(const toml::value &t, const char *table,
                                        const char *key)
{
   MFEM_VERIFY(t.contains(key), "[" << table << "] needs " << key << " (no default)");
   const toml::value &a = t.at(key);
   MFEM_VERIFY(a.is_array(), "[" << table << "] " << key << " must be an array of numbers");
   std::vector<real_t> out;
   for (const toml::value &v : a.as_array())
   {
      if (v.is_floating()) { out.push_back(static_cast<real_t>(v.as_floating())); }
      else if (v.is_integer()) { out.push_back(static_cast<real_t>(v.as_integer())); }
      else { MFEM_ABORT("[" << table << "] " << key << " must be an array of numbers"); }
   }
   return out;
}

/// An array of integers at key `key` of table `t`; required, may be empty.
inline std::vector<int> OptIntArray(const toml::value &t, const char *table, const char *key)
{
   MFEM_VERIFY(t.contains(key), "[" << table << "] needs " << key << " (no default)");
   const toml::value &a = t.at(key);
   MFEM_VERIFY(a.is_array(), "[" << table << "] " << key << " must be an array of integers");
   std::vector<int> out;
   for (const toml::value &v : a.as_array())
   {
      MFEM_VERIFY(v.is_integer(), "[" << table << "] " << key << " must be an array of integers");
      out.push_back(static_cast<int>(v.as_integer()));
   }
   return out;
}

} // namespace nitsche_opt
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_OPT_COMMON_HPP
