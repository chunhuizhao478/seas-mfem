// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qs_config_parser.hpp — TOML -> QsCaseConfig.  Tables:
//
//   [case]         name, source (the RACCOON input it reproduces)
//   [mesh]         type ("cartesian" | "slit"), nx, ny, xmin, xmax, ymin, ymax,
//                  slit_x, refine_steps,
//                  [[mesh.refine]]  lo, hi  |  center, length, width, height,
//                                   length_direction, width_direction
//                  [[mesh.nodeset]] name, boundary, lo, hi
//   [material]     E, nu
//   [phase_field]  model (AT1 | AT2), Gc, l, p, eta, decomposition,
//                  driving_average ("element" | "quadrature")
//   [[dirichlet]]  boundaries, component, value, rate
//   [reaction]     name, boundaries, component
//   [time]         dt, end_time, num_steps
//   [solver]       fixed_point_max_its, fixed_point_rel_tol,
//                  fixed_point_abs_tol, accept_on_max_fixed_point_iteration,
//                  automatic_scaling, nl_rel_tol, nl_abs_tol, linear_rel_tol,
//                  quadrature_order
//   [output]       dir, nodal_steps
//
// Unknown tables and keys throw.

#ifndef MFEM_PFF_QS_CONFIG_PARSER_HPP
#define MFEM_PFF_QS_CONFIG_PARSER_HPP

#include "qs_config.hpp"
#include "toml_reader.hpp"

#include <sstream>

namespace mfem
{
namespace pff
{

class QsCaseConfigParser : private TomlReader
{
public:
   static QsCaseConfig ParseFile(const std::string &path)
   {
      return Parse(toml::parse(path));
   }

   static QsCaseConfig ParseString(const std::string &s, const std::string &name = "string")
   {
      std::istringstream iss(s);
      return Parse(toml::parse(iss, name));
   }

private:
   static QsCaseConfig Parse(const toml::value &data)
   {
      QsCaseConfig c;
      Check(data, "root", {"case", "mesh", "material", "phase_field", "dirichlet", "reaction",
                           "time", "solver", "output"
                          });
      {
         const auto &t = data.at("case");
         Check(t, "case", {"name", "source"});
         c.name = Str(t, "name", "");
         c.source = Str(t, "source", "");
      }
      {
         const auto &t = data.at("mesh");
         Check(t, "mesh", {"type", "nx", "ny", "xmin", "xmax", "ymin", "ymax", "slit_x",
                           "refine_steps", "refine", "nodeset"
                          });
         c.mesh_type = Str(t, "type", "cartesian");
         if (c.mesh_type != "cartesian" && c.mesh_type != "slit")
         {
            throw std::invalid_argument("[mesh] type must be cartesian or slit");
         }
         c.nx = Int(t, "nx", 1);
         c.ny = Int(t, "ny", 1);
         c.xmin = Real(t, "xmin", 0.0);
         c.xmax = Real(t, "xmax", 1.0);
         c.ymin = Real(t, "ymin", 0.0);
         c.ymax = Real(t, "ymax", 1.0);
         c.slit_x = Real(t, "slit_x", c.xmin);
         c.refine_steps = Int(t, "refine_steps", 0);
         if (t.contains("refine"))
         {
            for (const auto &b : toml::get<toml::array>(t.at("refine")))
            {
               Check(b, "mesh.refine", {"lo", "hi", "center", "length", "width", "height",
                                        "length_direction", "width_direction"
                                       });
               RefineBox r;
               r.oriented = b.contains("center");
               if (r.oriented == (b.contains("lo") || b.contains("hi")))
               {
                  throw std::invalid_argument("[[mesh.refine]] give lo/hi or center/length/..., not both");
               }
               r.lo = Vec3(b, "lo", r.lo);
               r.hi = Vec3(b, "hi", r.hi);
               r.center = Vec3(b, "center", r.center);
               r.length = Real(b, "length", 0.0);
               r.width = Real(b, "width", 0.0);
               r.height = Real(b, "height", 0.0);
               r.length_direction = Vec3(b, "length_direction", r.length_direction);
               r.width_direction = Vec3(b, "width_direction", r.width_direction);
               c.refine.push_back(r);
            }
         }
         if (t.contains("nodeset"))
         {
            for (const auto &b : toml::get<toml::array>(t.at("nodeset")))
            {
               Check(b, "mesh.nodeset", {"name", "boundary", "lo", "hi"});
               NodeSetBox n;
               n.name = Str(b, "name", "");
               n.boundary = Str(b, "boundary", "");
               n.lo = Vec3(b, "lo", n.lo);
               n.hi = Vec3(b, "hi", n.hi);
               c.nodesets.push_back(n);
            }
         }
      }
      {
         const auto &t = data.at("material");
         Check(t, "material", {"E", "nu"});
         c.E = Real(t, "E", 0.0);
         c.nu = Real(t, "nu", 0.0);
      }
      {
         const auto &t = data.at("phase_field");
         Check(t, "phase_field", {"model", "Gc", "l", "p", "eta", "decomposition",
                                  "driving_average"
                                 });
         c.pf.crack.model = ParseCrackModel(Str(t, "model", "AT2"));
         c.pf.Gc = Real(t, "Gc", 0.0);
         c.pf.l = Real(t, "l", 0.0);
         c.pf.deg.p = Real(t, "p", 2.0);
         c.pf.deg.eta = Real(t, "eta", 1e-6);
         c.pf.driving = DrivingEnergy::Active;
         c.decomposition = ParseDecomposition(Str(t, "decomposition", "none"));
         const std::string avg = Str(t, "driving_average", "element");
         if (avg != "element" && avg != "quadrature")
         {
            throw std::invalid_argument("[phase_field] driving_average must be element or quadrature");
         }
         c.driving_element_average = (avg == "element");
      }
      if (data.contains("dirichlet"))
      {
         for (const auto &b : toml::get<toml::array>(data.at("dirichlet")))
         {
            Check(b, "dirichlet", {"boundaries", "component", "value", "rate"});
            DirichletSpec s;
            s.boundaries = Names(b, "boundaries");
            s.component = Int(b, "component", 0);
            s.value = Real(b, "value", 0.0);
            s.rate = Real(b, "rate", 0.0);
            c.dirichlet.push_back(s);
         }
      }
      if (data.contains("reaction"))
      {
         const auto &t = data.at("reaction");
         Check(t, "reaction", {"name", "boundaries", "component"});
         c.reaction_name = Str(t, "name", "F");
         c.reaction_boundaries = Names(t, "boundaries");
         c.reaction_component = Int(t, "component", 0);
      }
      {
         const auto &t = data.at("time");
         Check(t, "time", {"dt", "end_time", "num_steps"});
         c.dt = Real(t, "dt", 0.0);
         c.end_time = Real(t, "end_time", 0.0);
         c.num_steps = Int(t, "num_steps", 0);
         if (c.dt <= 0.0) { throw std::invalid_argument("[time] dt must be > 0"); }
      }
      if (data.contains("solver"))
      {
         const auto &t = data.at("solver");
         Check(t, "solver", {"fixed_point_max_its", "fixed_point_rel_tol", "fixed_point_abs_tol",
                             "accept_on_max_fixed_point_iteration", "automatic_scaling",
                             "nl_rel_tol", "nl_abs_tol",
                             "linear_rel_tol", "quadrature_order"
                            });
         c.fixed_point_max_its = Int(t, "fixed_point_max_its", c.fixed_point_max_its);
         c.fixed_point_rel_tol = Real(t, "fixed_point_rel_tol", c.fixed_point_rel_tol);
         c.fixed_point_abs_tol = Real(t, "fixed_point_abs_tol", c.fixed_point_abs_tol);
         c.accept_on_max_fixed_point_iteration =
            Bool(t, "accept_on_max_fixed_point_iteration", c.accept_on_max_fixed_point_iteration);
         c.automatic_scaling = Bool(t, "automatic_scaling", c.automatic_scaling);
         c.nl_rel_tol = Real(t, "nl_rel_tol", c.nl_rel_tol);
         c.nl_abs_tol = Real(t, "nl_abs_tol", c.nl_abs_tol);
         c.linear_rel_tol = Real(t, "linear_rel_tol", c.linear_rel_tol);
         c.quadrature_order = Int(t, "quadrature_order", c.quadrature_order);
      }
      if (data.contains("output"))
      {
         const auto &t = data.at("output");
         Check(t, "output", {"dir", "nodal_steps"});
         c.output_dir = Str(t, "dir", c.output_dir);
         if (t.contains("nodal_steps"))
         {
            for (const auto &v : toml::get<toml::array>(t.at("nodal_steps")))
            {
               c.nodal_output_steps.push_back(static_cast<int>(v.as_integer()));
            }
         }
      }
      return c;
   }
};

} // namespace pff
} // namespace mfem

#endif
