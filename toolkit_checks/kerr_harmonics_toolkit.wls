(* kerr_harmonics_toolkit.wls -- Black Hole Perturbation Toolkit spin-weighted spheroidal harmonics for the
   resonant modes of the photon on the Kerr light ring (Sec. IV of "Gravitational radiation from a photon on
   the light ring", E. Barausse).

   For a = 0 and for a = 1/2 and 9/10, prograde (sign = +1) and retrograde (sign = -1) light rings,
   l = 100, 200, 400, 800 and j = l - m = 0 ... 3, the mode (l, m) at the orbital frequency omega = m Omega = m/b (M = 1) has the
   spheroidicity gamma = a omega (negative for retrograde orbits, i.e. a -> -a with m > 0, the convention of
   the paper and of the notebook). The script stores the Toolkit's eigenvalue Lambda of
   SpinWeightedSpheroidalEigenvalue[-2, l, m, gamma] (the separation constant of the radial Teukolsky
   equation in the Toolkit's convention, the Lambda of Eq. (16) of the paper; Lambda = l(l+1) - s(s+1) at
   gamma = 0) and the equatorial values |S(Pi/2, 0)|^2 and |d_theta S(Pi/2, 0)|^2 of the harmonic
   SpinWeightedSpheroidalHarmonicS[-2, l, m, gamma], normalised to unity on the sphere (the function
   includes e^(i m phi); Int |S|^2 dOmega = 1). The spherical-expansion method is used (200 spherical
   harmonics; the eigenvalue and the two equatorial values are converged to about 10^-12 in relative terms,
   checked against 120 and 300 terms at l = 800, a = +-9/10; the Toolkit's default continued-fraction
   method does not converge at these large spheroidicities, gamma up to 253). At a = 0 the Toolkit's
   spin-weighted spherical harmonics SpinWeightedSphericalHarmonicY[-2, l, m, theta, phi] are used, with the
   eigenvalue l(l+1) - s(s+1) of SpinWeightedSpheroidalEigenvalue[-2, l, m, 0]: the paper's Lambda expansion
   and the equatorial limit reduce to the spherical results there (omega beta_b = m).

   These are the Toolkit data behind the two statements of Sec. IV "checked against the Toolkit's
   harmonics": Lambda = (m - a omega)^2 + (2j+1) omega beta_b + O(1) and |S(Pi/2)|^2 -> c_j
   Sqrt[omega beta_b]/(2 Pi^(3/2)), together with the odd-j analogue |S'(Pi/2)|^2 -> d_j (omega beta_b)^(3/2)
   /Pi^(3/2) and the ratios S'/S used by the source expansion of the notebook (Section IV). They are read
   by the notebook (Section IV, "Toolkit harmonics: the eigenvalue and the equatorial values") and by
   sympy_check/reproduce_paper.py.

   Run from this directory:   wolframscript -file kerr_harmonics_toolkit.wls      (about fifteen minutes)
   Needs the SpinWeightedSpheroidalHarmonics package of the Toolkit. Writes ../notebook/toolkit_harmonics.m
   (rows {a, sign, l, j, Lambda, S2, dS2} with a = |a| >= 0 and sign = +1 for a = 0, S2 = |S(Pi/2,0)|^2,
   dS2 = |d_theta S(Pi/2,0)|^2, machine precision, 80 rows) and prints every row. *)

<< SpinWeightedSpheroidalHarmonics`
scriptDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
outFile = FileNameJoin[{scriptDir, "..", "notebook", "toolkit_harmonics.m"}];

(* Equatorial circular null geodesics of Kerr (M = 1) with signed spin: aS < 0 is a retrograde orbit (a -> -a, m > 0). *)
rLR[aS_] := 2 (1 + Cos[2/3 ArcCos[-aS]]);
lrParams[aS_] := Module[{r0 = rLR[aS], D0, P0, b},
  D0 = r0^2 - 2 r0 + aS^2; P0 = 2 r0 D0/(r0 - 1);
  b = If[aS == 0, 3 Sqrt[3], (r0^2 + aS^2 - P0)/aS];
  <|"r0" -> r0, "b" -> b, "Omega" -> 1/b, "betab" -> Sqrt[b^2 - aS^2]|>];

nTerms = 200;
harmonics[aS_, l_, j_] := Module[{p = lrParams[aS], m = l - j, w, g, lam, S, S0, dS0},
  w = N[m p["Omega"], 30]; g = N[aS w];
  If[aS == 0,
    lam = SpinWeightedSpheroidalEigenvalue[-2, l, m, 0];
    S0 = SpinWeightedSphericalHarmonicY[-2, l, m, Pi/2, 0]; dS0 = Derivative[0, 0, 0, 1, 0][SpinWeightedSphericalHarmonicY][-2, l, m, Pi/2, 0];
    Return[<|"Lambda" -> N[lam], "S2" -> N[Abs[S0]^2], "dS2" -> N[Abs[dS0]^2], "gamma" -> 0., "wbetab" -> N[w p["betab"]]|>]];
  lam = Quiet[SpinWeightedSpheroidalEigenvalue[-2, l, m, g, Method -> {"SphericalExpansion", "NumTerms" -> nTerms}], {Eigensystem::arhm, Eigensystem::arh}];
  S = Quiet[SpinWeightedSpheroidalHarmonicS[-2, l, m, g, Method -> {"SphericalExpansion", "NumTerms" -> nTerms}], {Eigensystem::arhm, Eigensystem::arh}];
  S0 = S[Pi/2, 0]; dS0 = Derivative[1, 0][S][Pi/2, 0];
  <|"Lambda" -> lam, "S2" -> Abs[S0]^2, "dS2" -> Abs[dS0]^2, "gamma" -> g, "wbetab" -> N[w p["betab"]]|>];

rows = {};
Do[
  h = harmonics[sign a, l, j];
  Print["a = ", sign a, "  l = ", l, "  j = ", j, "  gamma = ", h["gamma"], "  Lambda = ", InputForm[h["Lambda"]],
    "  |S(Pi/2)|^2 = ", InputForm[h["S2"]], "  |S'(Pi/2)|^2 = ", InputForm[h["dS2"]],
    "   Lambda - (m - a w)^2 - (2j+1) w beta_b = ", h["Lambda"] - ((l - j) - h["gamma"])^2 - (2 j + 1) h["wbetab"]];
  AppendTo[rows, {a, sign, l, j, h["Lambda"], h["S2"], h["dS2"]}],
  {a, {0, 1/2, 9/10}}, {sign, If[a == 0, {1}, {1, -1}]}, {l, {100, 200, 400, 800}}, {j, 0, 3}];

(* write the table with a documenting comment on the first line (Get[] skips the comment) *)
str = OpenWrite[outFile];
WriteString[str, "(* toolkit_harmonics.m: output of toolkit_checks/kerr_harmonics_toolkit.wls. Rows {a, sign, l, j, Lambda, S2, dS2}: a = |a|/M (0, 1/2, 9/10), sign = +1 prograde / -1 retrograde (+1 at a = 0) (spheroidicity gamma = sign a omega, m = l - j > 0, omega = m/b), Lambda = SpinWeightedSpheroidalEigenvalue[-2, l, m, gamma] of the Black Hole Perturbation Toolkit (the Lambda of Eq. (16); l(l+1) - s(s+1) at gamma = 0), S2 = |S(Pi/2, 0)|^2 and dS2 = |d_theta S(Pi/2, 0)|^2 of SpinWeightedSpheroidalHarmonicS[-2, l, m, gamma] normalised to unity on the sphere (Int |S|^2 dOmega = 1, S including e^(i m phi)). M = 1, machine precision (spherical expansion with 200 terms; at a = 0 the Toolkit's SpinWeightedSphericalHarmonicY). *)\n"];
WriteString[str, ToString[rows, InputForm], "\n"];
Close[str];
Print["DONE: ", Length[rows], " rows written to ", outFile];
