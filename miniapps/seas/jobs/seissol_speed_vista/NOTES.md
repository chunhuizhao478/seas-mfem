# SeisSol decks for the Vista speed test: parity with the MFEM Nitsche runs

SeisSol source checked: `sst/ext/SeisSol` at a8e7f04f (paths below are
relative to its `src/`). MFEM paths are relative to `miniapps/seas/`.
The MFEM decks the speed jobs run are `mfem_decks/<prob>_nitsche_sst.toml`.
`check_parity.py` checks that each one parses to the same toml as its source
deck (`tpv205/configs/tpv205_nitsche.toml`, `tpv104/configs/tpv104_nitsche.toml`,
`tpv102/configs/tpv102_nitsche_pc_p2_kv0p1.toml`), except
`paraview_enabled`. The line numbers below are those of the `_sst` decks;
the source decks put the same keys 4 lines earlier.

Files:

| file | role |
|---|---|
| `decks/<prob>/parameters.par` | SeisSol namelist: FL, RS constants, t_0, mesh, CFL 0.5, ClusteredLTS 2, EndTime 12.0, output |
| `decks/<prob>/<prob>_fault.yaml` | easi: prestress, nuc_*, friction fields |
| `decks/<prob>/<prob>_material.yaml` | easi: rho, mu, lambda |
| `decks/<prob>/<prob>_faultreceivers.dat` | pickpoints = MFEM station CSV, same order |
| `decks/<prob>/<prob>_receivers.dat` | off-fault receivers = MFEM native surface list (empty for TPV205) |
| `check_parity.py` | MFEM deck + driver rules vs SeisSol deck + SeisSol rules, rel tol 1e-9 |
| `to_scec.py` | SeisSol `-faultreceiver-`/`-receiver-` files -> MFEM station file layout |
| `retag_msh2.py` | TPV205 only: boundary tags 101/103/105 -> 1/3/5 in a copy of the .msh, before pumgen |

Shared settings (all three decks): `OutputPointType = 3` (pickpoints only,
DRParameters.cpp:40-46, DRParameters.h:22-27); `SlipRateOutputType = 1`
(the default, DRParameters.cpp:67-70); `&Pickpoint printtimeinterval = 1`
(every time step of the cluster, OutputParameters.cpp:122), `OutputMask =
1 0 0 1 1 1 0 1 0 0 0 0` (SlipRate, FrictionAndState, TotalTractions, Slip,
AccumulatedSlip; order DataTypes.h:105-127, read at OutputParameters.cpp:126-128);
`ReceiverOutput = 1`, `pickdt = 0.005` (:145, :158); `Format = 10` (None,
OutputParameters.h:27); `WaveFieldOutput = 0`; `SurfaceOutput = 0`;
`EnergyOutput = 0` (DRParameters.cpp:154 for the DR part); `Checkpoint = 0`
(the default is on, OutputParameters.cpp:44); `OutputFile = 'output/<prob>'`
(the `output/` folder `run_pair_vista.sh` creates).
`XRef YRef ZRef = 0 -1e5 0`, `refPointMethod = 0` (DRParameters.cpp:32-38):
the plus side is y < 0, n = (0, 1, 0), strike = (1, 0, 0), dip = (0, 0, -1)
(MeshReader.cpp:139-163, Misc.cpp:19-34). This is MFEM's frame with
nhat = (0, -1, 0), shat = (1, 0, 0), dhat = (0, 0, -1)
(drivers/nitsche_dyn_driver.cpp:1647-1654), with the SeisSol normal = -nhat.

Material in all three: rho 2670, mu 32038120320, lambda 32043759360
(Cs 3464, Cp 6000). SeisSol names: ParameterDB.cpp:332-334.

## TPV205 (FL = 16, linear slip weakening)

| quantity | MFEM value | MFEM source | SeisSol value | SeisSol source |
|---|---|---|---|---|
| law | slip_weakening | `tpv205_nitsche_sst.toml:48` | FL = 16 | `parameters.par:35`, Misc.h:33-47 |
| lambda, mu, rho | 3.204375936e10, 3.203812032e10, 2670 | `:57`, `:58`, `:59` | same | `tpv205_material.yaml` |
| background stress | s_yy = -120e6, s_xy = 70e6 (tension +) | `:81` [stress] constant_tensor, driver:1567-1645 | same | `tpv205_fault.yaml` ConstantMap + LuaMap |
| patches | s_xy 81.6e6 (x in [-1.5, 1.5] km), 78.0e6 ([-9, -6] km), 62.0e6 ([6, 9] km), z in [-9, -6] km, closed boxes, last wins | [[stress_region]] x3 (`:88`, `:100`, `:112`), driver:1567-1645 | same, inclusive `>=`/`<=`, same order | `tpv205_fault.yaml` LuaMap |
| nucleation | none (patch 81.6 > 0.677 x 120 = 81.24 MPa) | no [nucleation] | nuc_* = 0 | `tpv205_fault.yaml`; required with nucleationcount 1, BaseDRInitializer.cpp:244-248, 374-441 |
| mu_s, mu_d, d_c, cohesion | 0.677, 0.525, 0.40 m, 0 | `:161`-`:164` | same | `tpv205_fault.yaml`; LinearSlipWeakeningInitializer.cpp:58-67 |
| forced rupture | none | - | none (no forced_rupture_time) | LinearSlipWeakening.h:198-207 |
| boundary tags | 103 fault, 101 free, 105 absorbing | `:167-169` | 3, 1, 5 after `retag_msh2.py` | FaceMap.cpp:98-117 |
| end time | 12 s | `:133` | 12.0 | `parameters.par:88` |
| mesh | `tpv205/mesh/tpv2053d_200m.geo` | - | `tpv2053d_200m.puml.h5` (pumgen of a retagged copy of the same .msh) | `parameters.par:64` |
| fault stations | 16 (incl. 7 at z = 0) | `tpv205/stations/tpv205_fault_stations.csv` | same order | `tpv205_faultreceivers.dat` |
| off-fault stations | none | dynamic/tpv205_setup.hpp:195-210 (intentionally empty) | empty file | `tpv205_receivers.dat` |
| KV damping | kv_eta 0.1 | `:128` | none | Open |
| Prakash-Clifton | v* 1, L 0.01 on sigma_n | `:193` | none | Open |
| CFL | 0.5 x C(p) | `:125` | CFL 0.5, ClusteredLTS 2 | `parameters.par:69-70` (`CFL`, `ClusteredLTS`) |

## TPV104 (FL = 103, slip law with strong rate weakening)

| quantity | MFEM value | MFEM source | SeisSol value | SeisSol source |
|---|---|---|---|---|
| law | rate_state / slip_law_strong_rate_weakening | `tpv104_nitsche_sst.toml:44`, `:124` | FL = 103 | `parameters.par:30` |
| lambda, mu, rho | as above | `:52`, `:53`, `:54` | same | `tpv104_material.yaml` |
| stress | s_yy = -120e6, s_xy = 40e6 | `:74` | same | `tpv104_fault.yaml` |
| f0, b, V0 | 0.6, 0.014, 1e-6 | `:121`, `:118`, `:120` | RS_f0, RS_b, RS_sr0 | `parameters.par:32-34`; DRParameters.cpp:104-108 |
| f_w | 0.2 | `:125` | RS_muW 0.2 | `parameters.par:35`; DRParameters.cpp:110-111 |
| V_init | 1e-16 | `:122` | RS_iniSlipRate1 1e-16, RS_iniSlipRate2 0 | `parameters.par:36-37` |
| L (Dc) | 0.4 m | `:119` | rs_sl0 0.4 | `tpv104_fault.yaml`; RateAndStateInitializer.cpp:139-151 |
| a | 0.01 inside, 0.02 outside, SCEC boxcar, x: 0 +- 15 km, z: -7.5 +- 7.5 km, transition 3 km | `:133` boxcar_taper; spatial_friction.cpp:2190-2205, spatial_friction.hpp:336-342 | rs_a LuaMap, same formula | `tpv104_fault.yaml` |
| V_w | 0.1 inside, 1.0 outside, same boxcar | `:135` | rs_srW LuaMap | `tpv104_fault.yaml`; RateAndStateInitializer.cpp:186 |
| psi_0 | InitialStatePsi(tau0, sigma0, V_init) | friction/dieterich_ruina.hpp:507-529 | FL = 103 init | RateAndStateInitializer.cpp:154-180 (equal at 2396 pts) |
| nucleation | dtau 45e6, compact bell R 3 km at (0, -7.5 km), T 1 s | `:155-156`; driver:1656-1684, spatial_nucleation.cpp:258-265, :20-29 | nuc_xy = 45e6 F(r), t_0 = 1.0 | `tpv104_fault.yaml`, `parameters.par:38`; GaussianNucleationFunction.h:22-31 |
| tags | 3 / 1 / 5 | `:159-161` | default face map | FaceMap.cpp:98-117 |
| end time | 12 s | `:87` | 12.0 | `parameters.par:88` |
| mesh | `tpv104/mesh/tpv104_200m.geo` | - | `tpv104_200m.puml.h5` | `parameters.par:64` |
| fault stations | 9 | `tpv104/stations/tpv104_fault_stations.csv` | same order | `tpv104_faultreceivers.dat` |
| off-fault stations | 6 (surf_0_9 ... surf_n12_n6), z = 0 | dynamic/tpv104_setup.hpp:209-219 | same order | `tpv104_receivers.dat` |
| KV / PC / CFL | 0.1 / v* 1, L 0.01 / 0.5 x C(p) | `:82` / `:185` / `:79` | none / none / 0.5 | Open |

## TPV102 (FL = 3, aging law)

| quantity | MFEM value | MFEM source | SeisSol value | SeisSol source |
|---|---|---|---|---|
| law | rate_state / aging_law | `tpv102_nitsche_sst.toml:42`, `:132` | FL = 3 | `parameters.par:31` |
| lambda, mu, rho | as above | `:51`, `:52`, `:53` | same | `tpv102_material.yaml` |
| stress | s_yy = -120e6, s_xy = 75e6 | `:73` | same | `tpv102_fault.yaml` |
| f0, b, V0 | 0.6, 0.012, 1e-6 | `:129`, `:126`, `:128` | RS_f0, RS_b, RS_sr0 | `parameters.par:33-35` |
| V_init | 1e-12 | `:130` | RS_iniSlipRate1 1e-12 | `parameters.par:36` |
| L (Dc) | 0.02 m | `:127` | rs_sl0 0.02 | `tpv102_fault.yaml` |
| a | 0.008 inside, 0.016 outside, same boxcar as TPV104 | `:137` | rs_a LuaMap | `tpv102_fault.yaml` |
| psi_0 / theta_0 | InitialStatePsi | dieterich_ruina.hpp:507-529 | FL = 3 theta_0 (psi = f0 + b ln(V0 theta / L)) | RateAndStateInitializer.cpp:104-132 (equal at 2396 pts) |
| nucleation | dtau 25e6, R 3 km, T 1 s | `:159-160` | nuc_xy = 25e6 F(r), t_0 = 1.0 | `tpv102_fault.yaml`, `parameters.par:38` |
| tags | 3 / 1 / 5 | `:163-165` | default face map | FaceMap.cpp:98-117 |
| end time | 12 s | `:87` | 12.0 | `parameters.par:88` |
| mesh | `tpv104/mesh/tpv104_200m.geo` (the gold run used it; tpv102_200m.geo differs only in comments) | - | `tpv104_200m.puml.h5` | `parameters.par:64` |
| stations | 9 fault, 6 surface | `tpv102/stations/tpv102_fault_stations.csv`, dynamic/tpv102_setup.hpp:178-189 | same order | `tpv102_faultreceivers.dat`, `tpv102_receivers.dat` |
| KV / PC / CFL | 0.1 / v* 1, L 0.01 / 0.5 x C(p) | `:82` / `:191` / `:78` | none / none / 0.5 | Open |

The TPV102 MFEM deck is the one of the p2 gold
(`tpv102/gold/nitsche-kv0p1_p2_cd_pc_200m_jobid7995488_git2e0c04a/`); its
physics equals `configs/tpv102_nitsche_pc_p2_kv0p1.toml` (only comments and
`[output.stations]` differ). There is no p3 gold for TPV102.

## Output mapping (to_scec.py)

SeisSol on-fault file `<OutputFile>-faultreceiver-NNNNN.dat`
(DynamicRupture/Output/OutputManager.cpp:320-337, header :346-452, rows
:556-571); NNNNN = line of the pickpoint file = CSV row. Columns are written
into MFEM's layout `time slip1 slip2 V1 V2 tau1 tau2 sigma_n <state>`
(dynamic/fault_station_file.hpp:320-327):

| MFEM column | SeisSol | note |
|---|---|---|
| slip1 (dip) | Sld | ReceiverBasedOutput.cpp:346-367 |
| slip2 (strike) | Sls | |
| V1, V2 | SRd, SRs | :438-447 |
| tau1, tau2 | Td0, Ts0 | total traction = initial + nucleation + change, :308-325 |
| sigma_n | -Pn0 | compression positive |
| 9th TPV205 `slip_path` | ASl | LSW has StV = 0 (Output/RateAndState.h:25-27) |
| 9th TPV104 `psi` | StV | |
| 9th TPV102 `log10_theta` | log10(StV) | StV is theta for FL = 3; MFEM writes log10((Dc/V0) exp((psi - f0)/b)), driver:2005-2017 |

Off-fault `<OutputFile>-receiver-NNNNN.dat` (ResultWriter/ReceiverWriter.cpp:93-149;
variables Equations/elastic/Model/Datastructures.h:36-37) -> `surface/<prob>_station_<name>.dat`,
`time vx vy vz` = Time, v1, v2, v3, the native MFEM surface layout
(dynamic/tpv104_setup.hpp:255-265). They are written to `surface/` so the
fault-station globs (`*_station_*.dat`) of the compare scripts do not take them.

## Assumptions

1. The meshes are the same Gmsh `.msh` files, converted by PUMGen
   (`build_speed_test_vista.sbatch`); the element set is not re-checked here.
2. LSW slip: MFEM's mu(D) uses |slip vector|; SeisSol uses the accumulated
   path (ASl). Equal while the slip direction does not turn; TPV205 is
   nearly unidirectional. Not compared.
3. Initial slip rate: SeisSol puts RS_iniSlipRate1 on the strike axis
   (tangent 1), MFEM along the initial shear. Here the shear is along strike,
   so they agree; magnitudes 1e-16 / 1e-12 m/s.
4. Cohesion: SeisSol strength = -cohesion - mu min(sigma_n, 0)
   (LinearSlipWeakening.h:158-160); 0 in both, so the sign does not matter.
5. TPV102 runs on the TPV104 200 m mesh, as the MFEM gold.
6. Pickpoints: SeisSol reports the value at the nearest fault Gauss point of
   the face that holds the point (OutputAux.cpp:209-227); MFEM uses the nearest
   node and averages copies. The 7 TPV205 stations at z = 0 were located by
   reproducing `isInsideFace` (PickPointBuilder.h:148-170, tolerance -1e-12)
   on a locally made `tpv2053d_200m.msh` (gmsh 4.12.1): all 16 found, with
   z = 0 giving exactly 0. A station not found aborts SeisSol (:184-209).
7. Receiver flush: `ReceiverOutputInterval` and the pickpoint `OutputInterval`
   stay at their defaults (1e100, OutputParameters.cpp:144, OutputParameters.h:23),
   so the data are held in memory and written at the end. This adds no
   LTS synchronisation points to the timed run.
8. The Lua bodies use branch form for the boxcar (|x| <= h -> 1, >= h + w -> 0,
   else the tanh formula), so the tanh term is never evaluated at its poles.
   check_parity.py evaluates the same Lua text.
9. Bell radius r is measured in the fault plane: r^2 = x^2 + (z + 7500)^2.
10. SeisSol's numerical flux and `GPwise` settings are left at their defaults.

## Open

1. **Kelvin-Voigt damping** (kv_eta 0.1 in all three MFEM decks): SeisSol
   a8e7f04f has none. The physics differ by this regularisation.
2. **Prakash-Clifton** (v* 1 m/s, L 0.01 m on sigma_n, all three decks):
   SeisSol has it only in FL = 6 (bimaterial LSW), where it acts on the
   strength (LinearSlipWeakening.h:280-295). This is a different law, so it
   is not used.
3. **CFL**: both decks say 0.5, but the definitions differ (MFEM 0.5 x C(p);
   SeisSol 0.5 with its own order-dependent factor and rate-2 LTS clusters).
   This changes time steps and speed, not the physics.
4. **TPV205 boundary tags**: `&MeshNml facemapfile` (MeshParameters.cpp:62-63)
   cannot be used at a8e7f04f: InitMesh.cpp:279-287 passes the file NAME to
   `YAML::Load`. Tags 101 and 105 would become dynamic-rupture faces
   (FaceMap.cpp:98-117), so `build_speed_test_vista.sbatch` renames the TPV205
   mesh's triangle tags 101/103/105 to 1/3/5 in a copy (`retag_msh2.py`) before
   pumgen, and records the resulting `tpv2053d_200m.puml.h5` in `MESHES.sha256`.
   check_parity.py checks the tags on that basis.
5. TPV205 has no MFEM off-fault stations, so the receivers file is empty and
   the deck sets `ReceiverOutput = 0` (an empty `RFileName` with receivers on
   is untested); `to_scec.py` makes no surface files for TPV205.
6. The MFEM Nitsche driver writes no off-fault stations. The 6 TPV104/TPV102
   receivers follow the native setups, so there is no MFEM Nitsche data to
   compare them with.
7. The fty namelist parser (SeisSol submodule) is not in the checked-out
   tree. Key spelling and case handling were checked against the `read*`
   calls in Initializer/Parameters/*.cpp, not against fty.
8. The SeisSol slip / slip-rate sign follows from the code (it matches
   MFEM's). It is not yet confirmed on run data. `to_scec.py` warns when the
   final h-slip and the initial h-shear have opposite signs.
9. `run_pair_vista.sh` runs `python3 -I check_parity.py`. `-I` ignores
   user site-packages, so PyYAML must be installed system-wide or in the
   active environment on Vista. Python must be 3.11 or later (tomllib).
