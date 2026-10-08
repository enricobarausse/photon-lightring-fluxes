(* schwarzschild_toolkit.wls -- Black Hole Perturbation Toolkit cross-checks for the photon on the
   Schwarzschild light ring.

   Paper:     Secs. II and III D of "Gravitational radiation from a photon on the light ring"
              (E. Barausse).

   Computes:  for l = 2, 3, 4, 5 and m = l, l - 1 the single-m fluxes to infinity and into the
              horizon (per E^2, M = 1) are computed in four ways and compared:
                RW-MST   photonMode of schwarzschild/photonLR.wl: the Zerilli/Regge-Wheeler fluxes
                         with the MST solutions of the Toolkit's ReggeWheeler package and the null
                         source of Sec. II;
                Teuk     the same null stress-energy tensor inserted in the Teukolsky equation
                         (s = -2) through the circular-orbit source of the Toolkit's Teukolsky
                         package, with E -> 1, L -> 3 Sqrt[3] and Upsilon_t -> r0^2 p^t = 27 (the
                         check quoted in Sec. II: agreement to the precision of the Toolkit's
                         solutions, nine to seventeen digits);
                NI       photonModeNI of schwarzschild/photonNI.wl, the arbitrary-precision direct
                         integration;
                MP       photonModeMP of schwarzschild/photonMP.wl, the machine-precision
                         ninth-order Runge-Kutta integrator used for the production runs
                         (agreement with the Toolkit to 10^-7 at l = 5 and 10^-5 at l = 2, Sec. III D).
   Inputs:    the ReggeWheeler and Teukolsky packages of the Toolkit; ../schwarzschild/photonNI.wl
              (which loads photonLR.wl) and ../schwarzschild/photonMP.wl.
   Outputs:   schwarzschild_toolkit_results.m next to this script (rows {l, m, method, Edot_I,
              Edot_H}, single m, 20 digits; the stored copy is read by sympy_check/reproduce_paper.py
              and by the notebook to verify the digits quoted in the paper); the fluxes and the
              relative differences with respect to RW-MST are also printed.
   Run:       from this directory:   wolframscript -file schwarzschild_toolkit.wls
   Cost:      a few minutes. *)

<< Teukolsky`
scriptDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
(* photonNI.wl loads photonLR.wl and the ReggeWheeler package *)
Get[FileNameJoin[{scriptDir, "..", "schwarzschild", "photonNI.wl"}]];
Get[FileNameJoin[{scriptDir, "..", "schwarzschild", "photonMP.wl"}]];

(* the photon orbit for the Toolkit's circular-orbit source: E = 1, L = b = 3 Sqrt[3], Upsilon_t = 27 *)
photonOrbit = <|"a" -> 0, "p" -> 3, "e" -> 0, "Inclination" -> 1, "Energy" -> 1,
   "AngularMomentum" -> 3 Sqrt[3],
   "Frequencies" -> <|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)" -> 27|>|>;
teukPhotonMode[l_, m_, wp_: 40] := Module[{w = N[m/(3 Sqrt[3]), wp], R, S, Z, assoc, fl},
  R = TeukolskyRadial[-2, l, m, 0, w, WorkingPrecision -> wp];
  S = SpinWeightedSpheroidalHarmonicS[-2, l, m, 0];
  Z = Teukolsky`ConvolveSource`Private`ConvolveSourcePointParticleCircular[-2, R, S, photonOrbit];
  assoc = <|"s" -> -2, "l" -> l, "m" -> m, "a" -> 0, "\[Omega]" -> w, "Eigenvalue" -> R["In"]["Eigenvalue"],
     "Amplitudes" -> Z, "Type" -> {"PointParticleCircular"}, "RadialFunctions" -> R, "AngularFunction" -> S,
     "Method" -> "MST"|>;
  fl = TeukolskyMode[assoc]["Fluxes"]["Energy"];
  <|"I" -> fl["\[ScriptCapitalI]"], "H" -> fl["\[ScriptCapitalH]"]|>];

(* relative difference, 3 digits *)
rel[x_, y_] := N[Abs[x/y - 1], 3];

res = {};
Do[
  mst = photonMode[l, m, WorkingPrecision -> 32]["Flux"];
  teuk = teukPhotonMode[l, m];
  ni = photonModeNI[l, m, 30]["Flux"];
  mp = photonModeMP[l, m]["Flux"];
  Print["l = ", l, " m = ", m];
  Print["  RW-MST  l*Edot_I = ", N[l mst["I"], 16], "  l*Edot_H = ", N[l mst["H"], 16]];
  Print["  Teuk    l*Edot_I = ", N[l teuk["I"], 16], "  l*Edot_H = ", N[l teuk["H"], 16],
    "   rel. diff. ", {rel[teuk["I"], mst["I"]], rel[teuk["H"], mst["H"]]}];
  Print["  NI      l*Edot_I = ", N[l ni["I"], 16], "  l*Edot_H = ", N[l ni["H"], 16],
    "   rel. diff. ", {rel[ni["I"], mst["I"]], rel[ni["H"], mst["H"]]}];
  Print["  MP      l*Edot_I = ", N[l mp["I"], 16], "  l*Edot_H = ", N[l mp["H"], 16],
    "   rel. diff. ", {rel[mp["I"], mst["I"]], rel[mp["H"], mst["H"]]}];
  AppendTo[res, {l, m, "RW-MST", N[mst["I"], 20], N[mst["H"], 20]}];
  AppendTo[res, {l, m, "Teuk", N[teuk["I"], 20], N[teuk["H"], 20]}];
  AppendTo[res, {l, m, "NI", N[ni["I"], 20], N[ni["H"], 20]}];
  AppendTo[res, {l, m, "MP", N[mp["I"], 20], N[mp["H"], 20]}],
  {l, {2, 3, 4, 5}}, {m, {l, l - 1}}];
Put[res, FileNameJoin[{scriptDir, "schwarzschild_toolkit_results.m"}]];
Print["DONE: ", Length[res], " rows written to schwarzschild_toolkit_results.m"];
