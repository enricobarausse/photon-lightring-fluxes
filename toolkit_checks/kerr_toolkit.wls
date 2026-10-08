(* kerr_toolkit.wls -- Black Hole Perturbation Toolkit reference fluxes for the photon on the Kerr
   light ring.

   Paper:     Sec. IV of "Gravitational radiation from a photon on the light ring" (E. Barausse).

   Computes:  the null stress-energy tensor of the photon is inserted in the Teukolsky equation
              (s = -2) through the circular-orbit source of the Toolkit's Teukolsky package (MST
              solutions, 40 digits), with E -> 1, L -> b (the impact parameter) and
              Upsilon_t -> r0^2 p^t (fluxes per E^2, M = 1), for a = 0, 1/2, 9/10 (prograde and
              retrograde), l = 2 ... 6 and j = l - m = 0, 1. These are the reference values of
              kerr_py/toolkit_reference.txt, against which kerr_py/validate_kerr.py checks the
              Python Teukolsky solver at l <= 5 (columns a, sign, l, j, l*Edot_I, l*Edot_H).
   Inputs:    the Teukolsky package of the Toolkit.
   Outputs:   kerr_toolkit_results.m next to this script (rows {a, sign, l, j, Edot_I, Edot_H},
              single m); every row is also printed.
   Run:       from this directory:   wolframscript -file kerr_toolkit.wls
   Cost:      about ten minutes. *)

<< Teukolsky`
scriptDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];

(* Equatorial circular null geodesics of Kerr (M = 1): sign = +1 prograde, -1 retrograde; E = 1, L = b. *)
rLR[a_, sign_] := 2 (1 + Cos[2/3 ArcCos[-sign a]]);
lrParams[a_, sign_] := Module[{r0 = rLR[a, sign], Dl, P, b, Om, Ut},
  Dl = r0^2 - 2 r0 + a^2;
  P = 2 r0 Dl/(r0 - 1);
  b = a + sign P/Sqrt[Dl];
  Om = 1/b;
  Ut = (r0^2 + a^2) P/Dl + a (b - a);
  <|"r0" -> r0, "b" -> b, "Omega" -> Om, "Ut" -> Ut, "P" -> P, "Delta" -> Dl|>];

Print["a = 0: ", N[lrParams[0, 1]], "  (r0 = 3, b = 3 Sqrt[3] = 5.196, Upsilon_t = 27)"];

(* check of the null geodesic: R(r0) = 0 and R'(r0) = 0 with R = ((r^2 + a^2) - a b)^2 - Delta (b - a)^2 *)
Do[p = lrParams[a, sg]; Rr = ((r^2 + a^2) - a p["b"])^2 - (r^2 - 2 r + a^2) (p["b"] - a)^2;
  Print["a = ", a, " sign = ", sg, "  {R, R'} at r0: ", N[{Rr, D[Rr, r]} /. r -> p["r0"], 20]],
  {a, {1/2, 9/10}}, {sg, {1, -1}}];

(* energy fluxes of the mode (l, m = sign (l - j)) from the Toolkit's circular-orbit source with the
   photon's orbit parameters (E = 1, L = b, Upsilon_t) *)
teukPhoton[a_, sign_, l_, j_, wp_: 40] := Module[{p = lrParams[a, sign], m, w, R, S, Z, orbit, assoc},
  m = sign (l - j); w = N[m p["Omega"], wp];
  orbit = <|"a" -> a, "p" -> N[p["r0"], wp], "e" -> 0, "Inclination" -> 1, "Energy" -> 1,
    "AngularMomentum" -> N[p["b"], wp],
    "Frequencies" -> <|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)" -> N[p["Ut"], wp]|>|>;
  R = TeukolskyRadial[-2, l, m, a, w, WorkingPrecision -> wp];
  S = SpinWeightedSpheroidalHarmonicS[-2, l, m, a w];
  Z = Teukolsky`ConvolveSource`Private`ConvolveSourcePointParticleCircular[-2, R, S, orbit];
  assoc = <|"s" -> -2, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> w, "Eigenvalue" -> R["In"]["Eigenvalue"],
     "Amplitudes" -> Z, "Type" -> {"PointParticleCircular"}, "RadialFunctions" -> R, "AngularFunction" -> S,
     "Method" -> "MST"|>;
  TeukolskyMode[assoc]["Fluxes"]["Energy"]];

res = {};
Do[
  fl = teukPhoton[a, sg, l, j];
  AppendTo[res, {a, sg, l, j, fl["\[ScriptCapitalI]"], fl["\[ScriptCapitalH]"]}];
  Print["a = ", a, " sign = ", sg, " l = ", l, " j = ", j,
    "  l*Edot_I = ", N[l fl["\[ScriptCapitalI]"], 12], "  l*Edot_H = ", N[l fl["\[ScriptCapitalH]"], 12]],
  {a, {0, 1/2, 9/10}}, {sg, If[a == 0, {1}, {1, -1}]}, {l, 2, 6}, {j, 0, 1}];
Put[res, FileNameJoin[{scriptDir, "kerr_toolkit_results.m"}]];
Print["DONE: ", Length[res], " rows written to kerr_toolkit_results.m"];
