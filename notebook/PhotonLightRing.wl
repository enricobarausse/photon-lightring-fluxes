(* ::Title:: *)
(*Gravitational radiation from a photon on the light ring*)

(* ::Subtitle:: *)
(*Companion notebook: multipole fluxes and the logarithmic divergence*)

(* ::Text:: *)
(*This notebook reproduces the calculations of the paper: the numbered equations, Eqs. (1)-(24), and the numbers quoted in the text, in Table I, in the figure captions and in the footnotes. Its sections follow those of the paper, and the text cell before each input cell says which equation or number of the paper it reproduces.*)

(* ::Subsubsection:: *)
(*Conventions*)

(* ::Text:: *)
(*Units are G = c = 1 and M = 1 for the mass of the black hole; the energy of the particle is E = 1, so that all fluxes are in units of E^2/M^2 (in Sections V and VI A, which deal with massive particles, the fluxes are per unit mu^2 where stated). Fluxes are for a single mode (l, m) with m > 0; the physical flux in the multipole l is twice the sum over m > 0, as stated after Eq. (4). The modes are labelled by j = l - m, so that j = 0 is the dominant mode l = m, even j the Zerilli (even-parity) sector and odd j the Regge-Wheeler (odd-parity) sector; a < 0 denotes retrograde orbits. The sources, master functions and amplitudes follow the conventions of the ReggeWheeler and Teukolsky packages of the Black Hole Perturbation Toolkit (bhptoolkit.org), as in the paper.*)

(* ::Subsubsection:: *)
(*Reading order*)

(* ::Text:: *)
(*The cells build on each other, so each section should be run after the previous ones: Section II defines the source and the jumps, III A-III C the large-l expansion and the coefficients kappa_j, III D the numerical integrator and the comparison with the stored run up to l = 12800, IV the Kerr calculation, V the comparison with the formulae of the 1970s, and VI A-VII the timelike orbits and the estimates of the Conclusions.*)

(* ::Subsubsection:: *)
(*Toolkit cells*)

(* ::Text:: *)
(*Five cells need the Black Hole Perturbation Toolkit (the ReggeWheeler and Teukolsky packages for the cross-checks at l <= 5, the SpinWeightedSpheroidalHarmonics package for the eigenvalue and harmonic checks and for the Kerr integrator). They are commented out, the text cell before each of them says so and states what it checks, and the rest of the notebook does not depend on them.*)

(* ::Subsubsection:: *)
(*Stored data*)

(* ::Text:: *)
(*The results of the long numerical runs are read from the repository: photon_big_results.m and timelike_results.m in the repository root, kerr_mathematica_runL1.m and kerr_mathematica_runL2.m in notebook/ (next to this notebook), and kerr_py/kerr_results.json. The repository root (repoDir, defined in the first cell) is taken to be the parent of the directory containing this notebook, or the current directory when the notebook is evaluated as a script.*)

(* ::Section:: *)
(*II. The problem and the formalism*)

(* ::Subsubsection:: *)
(*Eq. (1): the null geodesic*)

(* ::Text:: *)
(*The null geodesic of Eq. (1): the photon at r0 = 3 has p^t = E/f(r0) and p^phi = b E/r0^2; the null condition p.p = 0 fixes the impact parameter b = 3 Sqrt[3], and the orbital frequency is Omega = p^phi/p^t = 1/b. The stress-energy tensor of Eq. (1) has the form of that of a massive particle with mu u^alpha replaced by p^alpha, so the standard circular-orbit sources apply with E0 -> E and L0 -> b E.*)

ClearAll["Global`*"];
repoDir = Quiet[Check[ParentDirectory[NotebookDirectory[]], Directory[]]];   (* repository root: the parent of the notebook's directory, or the current directory *)
f[r_] := 1 - 2/r;
rstar[r_] := r + 2 Log[r/2 - 1];
bLR = 3 Sqrt[3]; r0LR = 3; OmegaLR = 1/bLR;
(* the null condition at r0 = 3 (E = 1): -f p^t p^t + r0^2 p^phi p^phi = 0 *)
Clear[b];
pt = 1/f[3]; pphi = b/3^2;
Print["b from p.p = 0: ", Select[b /. Solve[-f[3] pt^2 + 3^2 pphi^2 == 0, b], Positive], ",  Omega = p^phi/p^t = ", (pphi/pt) /. b -> bLR, " = 1/b"];
(* r0 = 3 is the circular null geodesic: the effective potential f/r^2 of null geodesics is stationary there *)
Print["d/dr (f/r^2) = 0 at r = ", r /. Solve[D[f[r]/r^2, r] == 0, r]];

(* ::Subsubsection:: *)
(*Eq. (2): source coefficients and jump conditions*)

(* ::Text:: *)
(*Eq. (2): the stress-energy coefficients A, C, F of the ReggeWheeler package of the Toolkit for a circular orbit, with E0 -> E and L0 -> b E for the photon, and the jump conditions that they produce on the master function Psi and on its r derivative at r0 (Sago, Nakano and Sasaki 2003; Toolkit conventions). The function jumpsZRW returns {[Psi], [dPsi/dr]}; Y and dY stand for the complex conjugates of Y_lm and of d_theta Y_lm at (theta, phi) = (Pi/2, 0). The same function serves a timelike circular orbit with E0 = E/mu and L0 = L/mu (fluxes then per unit mu^2), which is used in Sections V and VI A. The cell also defines the equatorial harmonics in a LogGamma form that remains accurate at large l, and checks it against the built-in harmonics.*)

(* jumps {[Psi], [dPsi/dr]} at r0 for a circular orbit; Y, dY = conj. of Y_lm and d_theta Y_lm at (Pi/2, 0) *)
jumpsZRW[l_, m_, r0_, E0_, L0_, Y_, dY_, even_: Automatic] := Module[{ev, EA, EC, EF, n, rm2M, np6M, term1, dterm1, coeff, term2},
  ev = If[even === Automatic, EvenQ[l + m], even];
  EA = -16 Pi (r0 - 2) E0/r0^3 Y;                                                    (* A of Eq. (2) *)
  EC = -16 Pi (r0 - 2) L0/r0^4/(l (l + 1)) dY;                                        (* C of Eq. (2) *)
  EF = -16 Pi (r0 - 2) L0^2/E0/r0^5/((l - 1) l (l + 1) (l + 2)) (l (l + 1) - 2 m^2) Y;  (* F of Eq. (2) *)
  If[ev,
    n = (l - 1) (l + 2); rm2M = r0 - 2; np6M = n r0 + 6;
    term1 = -r0^4 EA/2/np6M/rm2M;
    dterm1 = -r0^4 EA/2/np6M/rm2M (4/r0 - 1/rm2M - n/np6M);
    coeff = 2/r0/rm2M;
    term2 = r0^3/4 (r0^2 n (n - 2) + r0 (14 n - 36) + 96) EA/(rm2M np6M)^2 + (n + 2) r0^2/4 EF/rm2M;
    {term1, -coeff term1 - dterm1 + term2},
    {r0^3/(r0 - 2) EC, -2 r0^2/(r0 - 2)^2 EC + r0^2/(r0 - 2) EC - 3 r0^2/(r0 - 2) EC + r0^3/(r0 - 2)^2 EC}]];

(* exact harmonics for small l, and a LogGamma form of Y_{l,l-j}(Pi/2,0), d_theta Y_{l,l-j}(Pi/2,0) that is stable at large l *)
YLR[l_, m_] := SphericalHarmonicY[l, m, Pi/2, 0];
dYLR[l_, m_] := Derivative[0, 0, 1, 0][SphericalHarmonicY][l, m, Pi/2, 0];
YLRnum[l_, m_, prec_: 30] := Module[{j = l - m},
  If[EvenQ[j], N[(-1)^((2 l - j)/2) Sqrt[(2 l + 1)/(4 Pi)] Exp[(LogGamma[j + 1] + LogGamma[2 l - j + 1])/2 - l Log[2] - LogGamma[j/2 + 1] - LogGamma[(2 l - j)/2 + 1]], prec], 0]];
dYLRnum[l_, m_, prec_: 30] := Module[{j = l - m}, If[OddQ[j], N[Sqrt[j (2 l - j + 1)], prec] YLRnum[l, m + 1, prec], 0]];
(* check of the LogGamma form against the built-in harmonics *)
Print["max difference between the LogGamma form and the built-in harmonics, l <= 12: ", Max[Table[Abs[YLRnum[l, l - j] - YLR[l, l - j]] + Abs[dYLRnum[l, l - j] - dYLR[l, l - j]], {l, 2, 12}, {j, 0, l - 2}]]];

(* ::Subsubsection:: *)
(*Eqs. (3)-(4): master equation and flux normalization*)

(* ::Text:: *)
(*Eqs. (3)-(4): the master equation is solved in the tortoise coordinate with the Zerilli or Regge-Wheeler potential (defined in Section III B, where they are expanded) and the amplitudes Z are built from the two homogeneous solutions and the jumps as in Eq. (4). The flux normalization is the one stated after Eq. (4), Edot = c_pm omega^2 |Z|^2 with c_+ = (l-1)(l+2)/[4 Pi l(l+1)] and c_- = l(l+1)/[16 Pi (l-1)(l+2)]. In Eq. (4) the jumps are in the tortoise coordinate, [dPsi/dr*] = f(r0) [dPsi/dr], and the Wronskian W = 2 I omega A_inc makes the amplitudes independent of the normalization of the homogeneous solutions. The construction itself is carried out numerically in Section III D (photonMode).*)

cpm[l_, m_] := If[EvenQ[l + m], (l - 1) (l + 2)/(l (l + 1))/(4 Pi), l (l + 1)/((l - 1) (l + 2))/(16 Pi)];
fluxFromZ[l_, m_, w_, ZI_, ZH_] := <|"I" -> cpm[l, m] Abs[w ZI]^2, "H" -> cpm[l, m] Abs[w ZH]^2|>;

(* ::Section:: *)
(*III. Fluxes from a photon on the Schwarzschild light ring*)

(* ::Subsection:: *)
(*A. The source*)

(* ::Subsubsection:: *)
(*Eq. (5): the jumps at large l, and Eq. (6): the cancellation*)

(* ::Text:: *)
(*Eq. (5): the jumps for the photon (r0 = 3, b^2 = 27) with m = l - j, in the Zerilli sector (even j) and in the Regge-Wheeler sector (odd j); the odd-parity jumps are exact, and the cell checks all of them against the expressions of the paper. Eq. (6): for generic E0, L0 and r0 the O(l^0) part of the even-parity derivative jump is 4 Pi Y (L0^2/r0^3 - E0^2/(r0-2))/E0, proportional to p.p, and vanishes for the photon; for a massive particle it equals -4 Pi Y/(E0 r0), i.e. -4 Pi Y mu^2/(E r0). The cell also prints the two contributions separately, that of A (the energy density) and that of F (the azimuthal pressure), and the coefficients with which A and F enter the Zerilli jumps at leading order in l (text before Eq. 6).*)

Clear[l, m, j, Y, dY];
{PsiE, dPsiE} = jumpsZRW[l, m, 3, 1, 3 Sqrt[3], Y, dY, True] // Simplify;
Print["[Psi] even   = ", PsiE];
Print["[Psi'] even  = ", dPsiE];
Do[Print["  m = l - ", j, ":  [Psi'] = ", Series[dPsiE /. m -> l - j, {l, Infinity, 1}] // Normal // Simplify, "  (Eq. 5: -8 Pi (2j+1) Y/l)"], {j, 0, 2, 2}];
{PsiO, dPsiO} = jumpsZRW[l, m, 3, 1, 3 Sqrt[3], Y, dY, False] // Simplify;
Print["[Psi] odd    = ", PsiO, "   [Psi'] odd = ", dPsiO];
Print["differences from Eq. (5): even [Psi] - 8 Pi Y/(l(l+1)) = ", Simplify[PsiE - 8 Pi Y/(l (l + 1))], ",  odd [Psi] + 16 Sqrt[3] Pi dY/(l(l+1)) = ", Simplify[PsiO + 16 Sqrt[3] Pi dY/(l (l + 1))], ",  odd [Psi'] - 16 Pi dY/(Sqrt[3] l(l+1)) = ", Simplify[dPsiO - 16 Pi dY/(Sqrt[3] l (l + 1))]];
Print["tortoise-coordinate jumps: [dPsi/dr*] = f(3) [dPsi/dr] = [dPsi/dr]/", 1/f[3]];
(* Eq. (6): O(l^0) derivative jump for generic E0, L0 at r0 *)
Clear[E0, L0, r0];
dPsiGen = jumpsZRW[l, l, r0, E0, L0, Y, 0, True][[2]];
lead = Limit[dPsiGen, l -> Infinity, Assumptions -> r0 > 3] // Simplify;
Print["O(l^0) part of [Psi'] = ", lead, "  = 4 Pi Y (L0^2/r0^3 - E0^2/(r0-2))/E0: ", Simplify[lead - 4 Pi Y/E0 (L0^2/r0^3 - E0^2/(r0 - 2))] == 0, " (Eq. 6)"];
Print["massive particle (E0 = E/mu = (r0-2)/Sqrt[r0 (r0-3)], L0 = L/mu = r0/Sqrt[r0-3]): ", Simplify[lead /. {L0 -> Sqrt[r0^2/(r0 - 3)], E0 -> (r0 - 2)/Sqrt[r0 (r0 - 3)]}], "  = -4 Pi Y/(E0 r0), i.e. -4 Pi Y mu^2/(E r0) in units of mu"];
(* the two O(l^0) contributions separately: the A term (energy density, the E0^2 piece) and the F term (azimuthal pressure, the L0^2 piece) *)
leadA = lead /. L0 -> 0; leadF = Simplify[lead - leadA];
Print["A contribution: ", leadA, "   F contribution: ", leadF];
(* how A and F enter the Zerilli jumps (text before Eq. 6): coefficients of A and F in [Psi] and [Psi'], exact and at leading order in l *)
psiGen = jumpsZRW[l, l, r0, E0, L0, Y, 0, True][[1]];
EAgen = -16 Pi (r0 - 2) E0/r0^3 Y; EFgen = -16 Pi (r0 - 2) L0^2/E0/r0^5/((l - 1) l (l + 1) (l + 2)) (l (l + 1) - 2 l^2) Y;
coefA = Simplify[(dPsiGen /. L0 -> 0)/EAgen]; coefF = Simplify[(dPsiGen - (dPsiGen /. L0 -> 0))/EFgen];
Print["[Psi] = ", Simplify[psiGen/EAgen], " A  ~ ", Normal[Series[Simplify[psiGen/EAgen], {l, Infinity, 2}]], " A"];
Print["[Psi'] = (", coefA, ") A + (", coefF, ") F  ~ ", Normal[Series[coefA, {l, Infinity, 0}]], " A + ", Normal[Series[coefF, {l, Infinity, -2}]], " F"];
Print["for the photon (r0 = 3, E0 = 1, L0 = 3 Sqrt[3]): A term = ", leadA /. {r0 -> 3, E0 -> 1}, ",  F term = ", leadF /. {r0 -> 3, E0 -> 1, L0 -> 3 Sqrt[3]}, ",  sum = ", Simplify[(leadA + leadF) /. {r0 -> 3, E0 -> 1, L0 -> 3 Sqrt[3]}]];

(* ::Subsubsection:: *)
(*Equatorial harmonics at large l: c_j and d_j*)

(* ::Text:: *)
(*The equatorial harmonics at large l (Sec. III C): |Y_{l,l-j}(Pi/2,0)|^2 -> c_j Sqrt[l]/(2 Pi^(3/2)) with c_j = Binomial[j, j/2] 2^-j, and |d_theta Y_{l,l-j}|^2 -> d_j l^(3/2)/Pi^(3/2) with d_j = j Binomial[j-1, (j-1)/2] 2^(1-j). The closed form of |Y_{l,l-j}(Pi/2,0)|^2 is checked against the built-in harmonics and the limits are taken exactly.*)

cj[j_] := Binomial[j, j/2]/2^j;                     (* c_j, even j *)
dj[j_] := j Binomial[j - 1, (j - 1)/2]/2^(j - 1);   (* d_j, odd j *)
Ysq[l_, j_] := (2 l + 1)/(4 Pi) j! (2 l - j)!/(4^l ((j/2)!)^2 (((2 l - j)/2)!)^2);   (* |Y_{l,l-j}(Pi/2,0)|^2, j even *)
Print["closed form of |Y_{l,l-j}|^2 against the built-in harmonics, l <= 10: ", Max[Table[Abs[Ysq[l, j] - YLR[l, l - j]^2], {l, 2, 10}, {j, 0, l, 2}]]];
Print["|Y|^2 / (c_j Sqrt[l]/(2 Pi^(3/2))) -> ", Table[{j, Limit[Ysq[l, j]/(cj[j] Sqrt[l]/(2 Pi^(3/2))), l -> Infinity]}, {j, 0, 4, 2}]];
Print["|dY|^2 / (d_j l^(3/2)/Pi^(3/2)) -> ", Table[{j, Limit[(j (2 l - j + 1)) Ysq[l, j - 1]/(dj[j] l^(3/2)/Pi^(3/2)), l -> Infinity]}, {j, 1, 3, 2}]];

(* ::Subsection:: *)
(*B. The barrier top*)

(* ::Subsubsection:: *)
(*Eq. (7): expansion of the barrier about its top*)

(* ::Text:: *)
(*The Zerilli and Regge-Wheeler potentials of Eq. (3), expanded about r = 3 in the tortoise coordinate: V0, k and V3 at leading order in l, the distance epsilon = omega^2 - V0 of the emission frequency from the top, and the barrier parameter eta_j = -epsilon/Sqrt[2k] = j + 1/2 + O(1/l) of Eq. (7), for both parities. The cell also prints the two combinations used in Eqs. (8)-(9), omega/(2k)^(1/4) -> Sqrt[l/2] and omega (2k)^(1/4) -> Sqrt[2] l^(3/2)/27.*)

VRW[l_, r_] := f[r] (l (l + 1)/r^2 - 6/r^3);
VZ[l_, r_] := With[{lam = (l - 1) (l + 2)/2}, f[r] (2 lam^2 (lam + 1) r^3 + 6 lam^2 r^2 + 18 lam r + 18)/(r^3 (lam r + 3)^2)];
Vpot[l_, parity_, r_] := If[parity === "Zerilli", VZ[l, r], VRW[l, r]];
Dst[e_] := f[r] D[e, r];   (* d/dr* *)
Clear[l, j, r];
Do[
  V0 = Vpot[l, par, 3]; kk = -Dst[Dst[Vpot[l, par, r]]] /. r -> 3; V3 = Dst[Dst[Dst[Vpot[l, par, r]]]] /. r -> 3;
  Print[par, ": V0 = ", Series[V0, {l, Infinity, -1}] // Normal, ",  k = ", Series[kk, {l, Infinity, -1}] // Normal, ",  V3 = ", Series[V3, {l, Infinity, -1}] // Normal];
  eps = (l - j)^2/27 - V0;
  Print["   epsilon = ", Series[eps, {l, Infinity, 0}] // Normal // Simplify, ",   eta_j = -epsilon/Sqrt[2k] = ", Series[-eps/Sqrt[2 kk], {l, Infinity, 0}] // Normal // Simplify],
  {par, {"Zerilli", "ReggeWheeler"}}];
kZ = -Dst[Dst[VZ[l, r]]] /. r -> 3;
Print["omega/(2k)^(1/4) -> ", Series[((l - j)/Sqrt[27])/(2 kZ)^(1/4), {l, Infinity, 0}] // Normal // Simplify, ",   omega (2k)^(1/4) -> ", Series[((l - j)/Sqrt[27]) (2 kZ)^(1/4), {l, Infinity, -1}] // Normal // Simplify, "   (leading order)"];

(* ::Subsubsection:: *)
(*The cubic term and the parabolic region (footnote of Sec. III B)*)

(* ::Text:: *)
(*The cubic term and the parabolic region (footnote of Sec. III B and the paragraph that contains it). Writing omega^2 - V = p0^2 + delta with p0^2 = omega^2 - V0 + k x^2/2 and delta = -V3 x^3/6, the WKB phase Int Sqrt[p0^2 + delta] dx changes by Int delta/(2 p0) dx; for |x| >> k^(-1/4), where p0 -> Sqrt[k/2] x, the shift is -(V3/36) Sqrt[2/k] x^3. It is small for |x|^3 << Sqrt[k]/V3, which is proportional to 1/l, so the quadratic approximation is faithful in the parabolic region |x| << (Sqrt[k]/V3)^(1/3) ~ M l^(-1/3); the local condition |delta| << p0^2 only requires |x| << k/V3 ~ M. In the variable z = (2k)^(1/4) x of Weber's equation the cubic correction is V3 z^3/(6 (2k)^(5/4)) = O(l^(-1/2)) z^3 and the parabolic region is |z| << l^(1/6). The cell verifies these statements with the leading-order V0, k and V3 of the previous cell.*)

Clear[x, kk, V3];
(* first-order change of the WKB phase for |x| >> k^(-1/4), where p0 = Sqrt[k/2] x *)
phaseShift = Integrate[(-V3 x^3/6)/(2 Sqrt[kk/2] x), x];
Print["phase shift Int delta/(2 p0) dx = ", phaseShift, ";  difference from -(V3/36) Sqrt[2/k] x^3: ", Simplify[phaseShift + (V3/36) Sqrt[2/kk] x^3, kk > 0]];
(* leading-order coefficients of the barrier, from the Zerilli potential expanded in the previous cell *)
kLead = Limit[kZ/l^2, l -> Infinity] l^2; V3Lead = Limit[(Dst[Dst[Dst[VZ[l, r]]]] /. r -> 3)/l^2, l -> Infinity] l^2;
Print["k = ", kLead, ",  V3 = ", V3Lead, ",  Sqrt[k]/V3 = ", Simplify[Sqrt[kLead]/V3Lead, l > 0], " (proportional to 1/l):  parabolic region |x| << (Sqrt[k]/V3)^(1/3) = ", Simplify[(Sqrt[kLead]/V3Lead)^(1/3), l > 0], ";  local condition |x| << k/V3 = ", kLead/V3Lead];
Print["cubic term in z: V3/(2k)^(5/4) = ", Simplify[V3Lead/(2 kLead)^(5/4), l > 0], " = O(l^(-1/2));  the parabolic region in z: |z| << (2k)^(1/4) (Sqrt[k]/V3)^(1/3) = ", Simplify[(2 kLead)^(1/4) (Sqrt[kLead]/V3Lead)^(1/3), l > 0], " ~ l^(1/6)"];

(* ::Subsubsection:: *)
(*Distance from resonance: lambda_L = Omega*)

(* ::Text:: *)
(*Where the emission frequency sits (text after Eq. 7). The instability rate of the light ring is lambda_L = Sqrt[2k]/(2 omega) = Omega: for the parabolic barrier the quasinormal frequencies are the poles of the response, at eta = I (n + 1/2) (where Cosh[Pi eta] in |C|^2 vanishes), i.e. omega^2 - V0 = -I (n + 1/2) Sqrt[2k], so that Im omega_QNM = -(n + 1/2) Sqrt[2k]/(2 omega) and Re omega_QNM = Sqrt[V0] = (l + 1/2) Omega + O(1/l), the eikonal quasinormal frequency of Sec. II. The mode (l, l - j) therefore lies (j + 1/2) Omega = (j + 1/2) lambda_L below the real part of the fundamental mode, whatever the value of l, and for j = 0 this distance equals the damping rate of the fundamental mode.*)

Clear[n, wR];
Print["lambda_L = Sqrt[2k]/(2 omega) at leading order = ", Simplify[Sqrt[2 kLead]/(2 l/Sqrt[27]), l > 0], " = Omega = ", OmegaLR];
Print["Cosh[Pi eta] at eta = I (n + 1/2): ", Simplify[Cosh[Pi I (n + 1/2)], n \[Element] Integers]];
Print["Im omega_QNM from omega^2 - V0 = -I (n+1/2) Sqrt[2k] with omega = wR + I wI: ", wI /. Solve[2 wR wI == -(n + 1/2) Sqrt[2 kk], wI][[1]], ";  at leading order ", Simplify[-(n + 1/2) Sqrt[2 kLead]/(2 l/Sqrt[27]), l > 0], " = -(n + 1/2) Omega"];
Print["Re omega_QNM = Sqrt[V0] = ", Series[Sqrt[VZ[l, 3]], {l, Infinity, 0}] // Normal // Simplify, " + O(1/l) = (l + 1/2) Omega + O(1/l):  ", Simplify[Limit[(Sqrt[VZ[l, 3]] - (l + 1/2)/Sqrt[27]) l, l -> Infinity]], "/l"];
Print["distance of omega = (l - j) Omega from (l + 1/2) Omega: ", Simplify[((l + 1/2) - (l - j)) OmegaLR], " = (j + 1/2) lambda_L;  eta_j = (j + 1/2) in units of lambda_L"];

(* ::Subsubsection:: *)
(*Weber's equation and Eqs. (8)-(9)*)

(* ::Text:: *)
(*Weber's equation u'' + (z^2/4 - eta) u = 0 in z = (2k)^(1/4) x, and its purely transmitted solution u = C D_nu(-E^(-I Pi/4) z), nu = -1/2 - I eta. With D_nu(0) = 2^(nu/2) Sqrt[Pi]/Gamma[(1-nu)/2], D_nu'(0) = -2^((nu+1)/2) Sqrt[Pi]/Gamma[-nu/2] and |C|^2 = omega E^(-Pi eta/2)/((2k)^(1/4) Cosh[Pi eta]), the values at the top are Eqs. (8)-(9); d/dr* = (2k)^(1/4) d/dz. This cell defines F(eta) and the reduced quantities u0sqA = |u(0)|^2/(omega/(2k)^(1/4)) and Du0sqA = |d_r* u(0)|^2/(omega (2k)^(1/4)) used throughout, checks that D_nu(-E^(-I Pi/4) z) solves Weber's equation, and checks Eqs. (8)-(9) against D_nu(0), D_nu'(0) and |C|^2.*)

G34[x_] := Abs[Gamma[3/4 + I x/2]]^2; G14[x_] := Abs[Gamma[1/4 + I x/2]]^2;
Fpar[x_] := Exp[-Pi x/2]/(Cosh[Pi x] G34[x]);                     (* F(eta) of the paper *)
u0sqA[x_] := Pi Fpar[x]/Sqrt[2];                                    (* Eq. (8): |u(0)|^2 / (omega/(2k)^(1/4)) *)
Du0sqA[x_] := Sqrt[2] Pi Exp[-Pi x/2]/(Cosh[Pi x] G14[x]);          (* Eq. (9): |u'(0)|^2 / (omega (2k)^(1/4)) *)
(* the transmitted solution and Weber's equation *)
Clear[z, eta, w, k];
uD[eta_, z_] := ParabolicCylinderD[-1/2 - I eta, -Exp[-I Pi/4] z];
Print["Weber residual of D_nu(-E^(-I Pi/4) z) at eta = 1/2, z = 1.3: ", Chop[N[(D[uD[eta, z], {z, 2}] + (z^2/4 - eta) uD[eta, z]) /. {eta -> 1/2, z -> 13/10}, 20]]];
(* Eqs. (8)-(9) from D_nu(0), D_nu'(0) and |C|^2 *)
nu = -1/2 - I eta;
Dnu0 = 2^(nu/2) Sqrt[Pi]/Gamma[(1 - nu)/2]; dDnu0 = -2^((nu + 1)/2) Sqrt[Pi]/Gamma[-nu/2];
Print["D_nu(0), D_nu'(0) against the built-in function at eta = 1/2: ", Chop[N[{Dnu0 - ParabolicCylinderD[nu, 0], dDnu0 - Derivative[0, 1][ParabolicCylinderD][nu, 0]} /. eta -> 1/2, 20]]];
C2 = w Exp[-Pi eta/2]/((2 k)^(1/4) Cosh[Pi eta]);
u0sq = C2 Abs[Dnu0]^2; Du0sq = C2 Sqrt[2 k] Abs[dDnu0]^2;      (* |u(0)|^2 and |d_r* u(0)|^2 *)
Print["Eq. (8) and Eq. (9), ratio of |C|^2 |D_nu(0)|^2 and |C|^2 Sqrt[2k] |D_nu'(0)|^2 to the right-hand sides, eta = 1/2, 3/2, 5/2 (w = k = 1): ",
  Table[N[{u0sq/((w/(2 k)^(1/4)) u0sqA[eta]), Du0sq/(w (2 k)^(1/4) Du0sqA[eta])} /. {eta -> x, w -> 1, k -> 1}, 15], {x, {1/2, 3/2, 5/2}}]];

(* ::Subsubsection:: *)
(*Normalization |C|^2 and the transmission coefficient*)

(* ::Text:: *)
(*Numerical check of the normalization |C|^2 and of the transmission coefficient |T|^2 = 1/(1 + E^(2 Pi eta)), which is 0.041 at eta = 1/2 (the barrier reflects 96% of the m = l waves). The solution is decomposed at z = 40 into incident and reflected WKB waves (iterated Riccati form of the local wavenumber); |u(0)|^2 for unit incident amplitude, in units of omega/(2k)^(1/4), must equal u0sqA. This is done first with the exact parabolic-cylinder function and then with the model equation u'' + (z^2/4 - eta - g z^3) u = 0 integrated numerically (u0sqModel), which is used in Section III C with g different from zero for the cubic correction of the barrier. The agreement is at the level of 10^-6, the accuracy of the WKB boundary data at |z| = 40.*)

(* local WKB wavenumber, iterated Riccati form, for u'' + (z^2/4 - eta - g z^3) u = 0 *)
QWeber[eta_, g_, sign_] := Module[{z, p, q}, p = Sqrt[z^2/4 - eta - g z^3]; q = sign p; Do[q = sign Sqrt[p^2 + I D[q, z]], {3}]; Function @@ {z, q}];
(* incident amplitude of a solution with values (uB, duB) at z = Z, for a wave incident from the right ("In") or from the left ("Up") *)
incAmp[eta_, g_, side_, Z_, uB_, duB_] := Module[{Qp = QWeber[eta, g, +1][Z], Qm = QWeber[eta, g, -1][Z]},
  If[side == "In", (duB - I Qp uB)/(I (Qm - Qp)), (duB - I Qm uB)/(I (Qp - Qm))]];
(* exact solution: |u(0)|^2 and |T|^2 for unit incident WKB amplitude from z = +Z, transmitted toward z = -Z *)
u0sqExact[eta_, Z_: 40] := Module[{p, Ainc, pB, pA, duD},
  duD[zz_] := -Exp[-I Pi/4] Derivative[0, 1][ParabolicCylinderD][-1/2 - I eta, -Exp[-I Pi/4] zz];
  p[zz_] := Sqrt[zz^2/4 - eta]; pB = p[Z]; pA = p[-Z];
  Ainc = incAmp[eta, 0, "In", Z, uD[eta, N[Z, 30]], duD[N[Z, 30]]];
  <|"u0sq" -> Abs[uD[eta, 0]]^2/(Abs[Ainc]^2 pB), "T2" -> Abs[uD[eta, N[-Z, 30]]]^2 pA/(Abs[Ainc]^2 pB)|>];
Print["exact D_nu, eta = 1/2: |u(0)|^2 / u0sqA = ", u0sqExact[1/2]["u0sq"]/u0sqA[1/2], ",  |T|^2 = ", u0sqExact[1/2]["T2"], " against 1/(1 + E^Pi) = ", N[1/(1 + Exp[Pi])]];
Print["exact D_nu, eta = 3/2: |u(0)|^2 / u0sqA = ", u0sqExact[3/2]["u0sq"]/u0sqA[3/2], ",  |T|^2 = ", u0sqExact[3/2]["T2"], " against ", N[1/(1 + Exp[3 Pi])]];
(* numerical model: u'' + (z^2/4 - eta - g z^3) u = 0; side "In": unit incident from z -> +Infinity, transmitted to the left; "Up": the mirror image *)
u0sqModel[eta_, g_, side_, Z_: 40] := Module[{z, zz, zA, zB, dir, Q, sol, Psi, Ainc, pB, p},
  p = Sqrt[z^2/4 - eta - g z^3];
  If[side == "In", zA = -Z; zB = Z; dir = -1, zA = Z; zB = -Z; dir = +1];
  Q = QWeber[eta, g, dir];
  sol = NDSolveValue[{Psi''[zz] + (zz^2/4 - eta - g zz^3) Psi[zz] == 0, Psi[zA] == 1, Psi'[zA] == I Q[zA]}, Psi, {zz, Min[zA, zB], Max[zA, zB]},
     PrecisionGoal -> 12, AccuracyGoal -> 12, MaxSteps -> 10^6, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9, "StiffnessTest" -> False}];
  Ainc = incAmp[eta, g, side, zB, sol[zB], sol'[zB]]; pB = p /. z -> zB;
  (* WKB flux of the incident wave is pB |Ainc|^2, of the transmitted wave p(zA) (unit amplitude there) *)
  <|"u0sq" -> Abs[sol[0]]^2/(Abs[Ainc]^2 pB), "Du0sq" -> Abs[sol'[0]]^2/(Abs[Ainc]^2 pB), "T2" -> (p /. z -> zA)/(Abs[Ainc]^2 pB)|>];
Print["model equation, eta = 1/2: |u(0)|^2 / u0sqA = ", u0sqModel[1/2, 0, "In"]["u0sq"]/u0sqA[1/2], ",  |u'(0)|^2 / Du0sqA = ", u0sqModel[1/2, 0, "In"]["Du0sq"]/Du0sqA[1/2], ",  |T|^2 = ", u0sqModel[1/2, 0, "In"]["T2"]];
Print["model equation, eta = 3/2: |u(0)|^2 / u0sqA = ", u0sqModel[3/2, 0, "In"]["u0sq"]/u0sqA[3/2]];

(* ::Subsection:: *)
(*C. The fluxes*)

(* ::Subsubsection:: *)
(*Eqs. (10)-(13): kappa_j and kappa*)

(* ::Text:: *)
(*Eqs. (10)-(13). The flux formula after Eq. (9), Edot = (c_pm/4) |u(0) [dPsi/dr*] -+ [Psi] d_r* u(0)|^2, with the leading terms: for even j, kappa_j = l |u(0)|^2 [dPsi/dr*]^2 c_+/4 with [dPsi/dr*] = -8 Pi (2j+1) Y/(3 l), omega/(2k)^(1/4) -> Sqrt[l/2] and |Y|^2 -> c_j Sqrt[l]/(2 Pi^(3/2)); for odd j, kappa_j = l |d_r* u(0)|^2 [Psi]^2 c_-/4 with [Psi] = -16 Sqrt[3] Pi dY/l^2, omega (2k)^(1/4) -> Sqrt[2] l^(3/2)/27 and |dY|^2 -> d_j l^(3/2)/Pi^(3/2). The results are compared with the closed forms of Eqs. (11)-(12), and the values kappa_0 = 0.0277668, kappa_1 = 0.0037670, kappa_2 = 2.7521 10^-4, kappa_3 = 1.6604 10^-5, kappa_4 = 9.23 10^-7, Sum_j kappa_j = 0.0318265 and kappa = 2 Sum_j kappa_j = 0.063653 of Eq. (13) are obtained.*)

kappaEven[j_] := Simplify[(1/(16 Pi)) u0sqA[j + 1/2] Sqrt[1/2] (64 Pi^2 (2 j + 1)^2/9) cj[j]/(2 Pi^(3/2))];
kappaOdd[j_] := Simplify[(1/(64 Pi)) Du0sqA[j + 1/2] (Sqrt[2]/27) (768 Pi^2) dj[j]/Pi^(3/2)];
kappaJ[j_] := If[EvenQ[j], kappaEven[j], kappaOdd[j]];
Print["kappa_0 = ", kappaEven[0], " = ", N[kappaEven[0], 8]];
Print["kappa_1 = ", kappaOdd[1], " = ", N[kappaOdd[1], 8]];
Print["differences from the closed forms of Eqs. (11)-(12), j = 0..5: ",
  Table[N[kappaJ[j] - If[EvenQ[j], Sqrt[Pi]/9 (2 j + 1)^2 cj[j] Fpar[j + 1/2], 8 Sqrt[Pi]/9 dj[j] Exp[-Pi (j + 1/2)/2]/(Cosh[Pi (j + 1/2)] G14[j + 1/2])], 20], {j, 0, 5}]];
TableForm[Table[{j, N[kappaJ[j], 8], N[kappaJ[j]/kappaJ[0], 6]}, {j, 0, 5}], TableHeadings -> {None, {"j", "kappa_j", "kappa_j/kappa_0"}}]
sumKappa = N[Sum[kappaJ[j], {j, 0, 20}], 12];
kappaTot = 2 sumKappa;
Print["Sum_j kappa_j = ", sumKappa, "    kappa = 2 Sum_j kappa_j = ", kappaTot, "   (Eq. 13: 0.0318265, 0.063653)"];
Print["kappa_{j+1}/kappa_j for j = 0..4: ", N[Table[kappaJ[j + 1]/kappaJ[j], {j, 0, 4}], 3], "  (the kappa_j fall off roughly like E^(-Pi j), E^-Pi = ", N[Exp[-Pi], 3], ")"];
Print["fraction of Sum_j kappa_j in the terms j >= 5 (not captured by j <= 4, Sec. III D): ", N[Sum[kappaJ[j], {j, 5, 20}]/sumKappa, 3], ";  in the terms j >= 4 (omitted in the Kerr comparison, Sec. IV): ", N[Sum[kappaJ[j], {j, 4, 20}]/sumKappa, 3]];

(* ::Subsubsection:: *)
(*Eq. (14): the phase ratio*)

(* ::Text:: *)
(*Eq. (14): the ratio u(0)/d_r* u(0) at the barrier top, (2k)^(1/4) times E^(I Pi/4) Gamma[1/4 + I eta/2]/(Sqrt[2] Gamma[3/4 + I eta/2]), whose phase sets the cross term; its modulus is consistent with Eqs. (8)-(9).*)

phaseRatio[eta_] := Exp[I Pi/4] Gamma[1/4 + I eta/2]/(Sqrt[2] Gamma[3/4 + I eta/2]);   (* (2k)^(1/4) u(0)/u'(0) *)
Print["(2k)^(1/4) u(0)/u'(0) at eta = 1/2: ", N[phaseRatio[1/2], 12], ",  phase = ", N[Arg[phaseRatio[1/2]], 10]];
(* consistency with Eqs. (8)-(9): |u(0)|^2 / |u'(0)|^2 = |phaseRatio|^2 / Sqrt[2k], i.e. u0sqA/Du0sqA = |phaseRatio|^2 *)
Print["u0sqA/Du0sqA - |phaseRatio|^2 at eta = 1/2: ", Chop[N[u0sqA[1/2]/Du0sqA[1/2] - Abs[phaseRatio[1/2]]^2]]];

(* ::Subsubsection:: *)
(*The O(l^(-1/2)) asymmetry sigma_j*)

(* ::Text:: *)
(*The O(l^(-1/2)) asymmetry sigma_j between the two fluxes, Eq. (10) and the remarks after Eq. (13). sigma_j has two contributions: the cross term of the flux formula, fixed by the phase of Eq. (14), and the cubic term of the barrier, which makes |v(0)| differ from |u(0)| (for odd j, |d_r* v(0)| from |d_r* u(0)|). The cubic contribution is obtained from the model equation u'' + (z^2/4 - eta - g z^3) u = 0 with g = V3/(6 (2k)^(5/4)) = 0.068 l^(-1/2), by integrating from both sides at a few small g and extrapolating to g -> 0 (g must stay well below 1/(4 Z) = 0.006 for the WKB data at |z| = 40 to be valid). For j = 0 the cross term is 1.039 and the cubic term -0.231, so that sigma_0 = 0.81; the sigma_j with j >= 1 are negative. Weighting the sigma_j with the kappa_j gives the m-summed asymmetry (Edot_I - Edot_H)/(Edot_I + Edot_H) -> sigma_bar/Sqrt[l] with sigma_bar = 0.61 (displayed after Eq. 14).*)

(* cross term, relative to the leading term kappa_j: even j from |u(0)|^2 [Psi']^2, odd j from |u'(0)|^2 [Psi]^2; the factor 1/Sqrt[27] is 1/(2k)^(1/4) Sqrt[l] at leading order *)
crossTerm[j_] := Module[{x = j + 1/2, ph},
  ph = Arg[Exp[I Pi/4] Gamma[1/4 + I x/2]/Gamma[3/4 + I x/2]];
  If[EvenQ[j],
    N[(-2 Sqrt[u0sqA[x] Du0sqA[x]] (1/Sqrt[27]) Cos[ph] (-8 Pi (2 j + 1)/3) (8 Pi) cj[j]/(2 Pi^(3/2))/(16 Pi))/kappaEven[j], 10],
    N[(-2 Sqrt[u0sqA[x] Du0sqA[x]] (1/Sqrt[27]) Cos[ph] (-256 Pi^2/3) dj[j]/Pi^(3/2)/(64 Pi))/kappaOdd[j], 10]]];
Print["cross term, j = 0: ", crossTerm[0], "  (paper: 1.039);  j = 1: ", crossTerm[1], ";  j = 2: ", crossTerm[2]];
(* cubic term of the barrier: V = l^2 g0(r) + O(l), g0 = f/r^2, so g = V3/(6 (2k)^(5/4)) is proportional to l^(-1/2) *)
g0[r_] := f[r]/r^2;
V2c = Dst[Dst[g0[r]]] /. r -> 3; V3c = Dst[Dst[Dst[g0[r]]]] /. r -> 3;
gcoef = V3c/(6 (2 (-V2c))^(5/4));
Print["g = ", N[gcoef], " l^(-1/2)"];
gList = {0.001, 0.002, 0.003, 0.005};
cubicTerm[j_] := Module[{x = j + 1/2, key = If[EvenQ[j], "u0sq", "Du0sq"], asym},
  asym = Table[With[{uIn = u0sqModel[x, g, "In"][key], uUp = u0sqModel[x, g, "Up"][key]}, {g, (uIn - uUp)/(uIn + uUp)/g}], {g, gList}];
  (Fit[asym, {1, g, g^2}, g] /. g -> 0) N[gcoef]];
asym0 = Table[With[{uIn = u0sqModel[1/2, g, "In"]["u0sq"], uUp = u0sqModel[1/2, g, "Up"]["u0sq"]}, {g, (uIn - uUp)/(uIn + uUp)/g}], {g, gList}];
Print["j = 0: (|u(0)|^2 - |v(0)|^2)/(|u(0)|^2 + |v(0)|^2)/g at g = ", gList, ": ", asym0[[All, 2]]];
cubic0 = cubicTerm[0];
Print["cubic contribution to sigma_0 = ", cubic0, "  (paper: -0.231);  sigma_0 = ", crossTerm[0] + cubic0, "  (paper: 0.81)"];
sigmaJ = Table[{j, crossTerm[j], cubicTerm[j]}, {j, 0, 4}];
TableForm[Map[{#[[1]], #[[2]], #[[3]], #[[2]] + #[[3]]} &, sigmaJ], TableHeadings -> {None, {"j", "cross term", "cubic term", "sigma_j"}}]
sigBar = Sum[(sigmaJ[[j + 1, 2]] + sigmaJ[[j + 1, 3]]) N[kappaJ[j]], {j, 0, 4}]/N[Sum[kappaJ[j], {j, 0, 4}]];
Print["m-summed asymmetry sigma_bar = Sum_j kappa_j sigma_j / Sum_j kappa_j = ", sigBar, "  (paper: 0.61)"];
Print["relative difference of the two fluxes at l = 100, 2 sigma_bar/Sqrt[l] = ", N[2 sigBar/10, 3], "  (paper: 12%)"];

(* ::Subsubsection:: *)
(*The horizon fraction*)

(* ::Text:: *)
(*The horizon fraction (end of Sec. III C): with the asymptotic formulae Edot^{inf,H} = kappa/l (1 +- sigma_bar/Sqrt[l]) summed from l = 2 to l_max, the fraction of the total that falls into the horizon approaches one half only logarithmically, because the l^(-1/2) terms cancel in the sum of the two fluxes but not in their ratio: it is 40% for l_max = 100-140 (Gundlach et al. 2012, truncating at l = 140 and extrapolating in radius to the light ring, found 42%).*)

horizonFraction[lmax_] := With[{l = Range[2, lmax]}, Total[(1 - sigBar/Sqrt[l])/l]/Total[2/l]];
TableForm[Table[{lmax, horizonFraction[lmax]}, {lmax, {20, 50, 100, 140, 400, 10^4, 10^6}}], TableHeadings -> {None, {"l_max", "horizon fraction"}}]

(* ::Subsection:: *)
(*D. Numerical check (Fig. 1)*)

(* ::Subsubsection:: *)
(*The integrator*)

(* ::Text:: *)
(*The Green's function construction of Sec. II carried out numerically, without any large-l approximation. For each (l, m) the homogeneous Zerilli or Regge-Wheeler equation is integrated in r* with a ninth-order explicit Runge-Kutta method, from r* = -40 with purely ingoing data (Psi_in) and from r = 40 with purely outgoing data (Psi_up); the boundary data are the iterated Riccati-WKB waves, with the amplitude integrals that make them unit incident waves. photonMode returns the fluxes to infinity and into the horizon from Eq. (4) with the jumps of Section II, and the unit-incident values u(0), d_r* u(0), v(0), d_r* v(0) at r0 = 3. Machine precision.*)

QWKB[l_, parity_, w_, sign_, niter_: 3] := Module[{r, p, Q},
  p = Sqrt[w^2 - Vpot[l, parity, r]]; Q = sign p;
  Do[Q = sign Sqrt[p^2 + I f[r] D[Q, r]], {niter}];
  Function @@ {r, Q}];

photonMode[l_, m_, rsA_: -40, rB_: 40, pg_: 12] := Module[{parity, w, Qin, Qup, rA, psiA, dpsiA, psiB, dpsiB, solIn, solUp, V, eqs, ampA, ampB, rsB, r, Psi, rr, PsiIn, dPsiIn, PsiUp, dPsiUp, W, J, P, dP, ZI, ZH, Ainc},
  parity = If[EvenQ[l + m], "Zerilli", "ReggeWheeler"];
  w = N[m OmegaLR];
  V = Vpot[l, parity, #] &;
  rA = N[r /. FindRoot[rstar[r] == rsA, {r, 2 + 2 Exp[(rsA - 2)/2], 2 + 10^-30, 3}, WorkingPrecision -> 40, AccuracyGoal -> 35], 40];
  Qin = QWKB[l, parity, w, -1]; Qup = QWKB[l, parity, w, +1];
  ampA = Exp[-NIntegrate[Im[Qin[rr]]/f[rr], {rr, 2, rA}, PrecisionGoal -> 12, MaxRecursion -> 30]];
  ampB = Exp[NIntegrate[Im[Qup[rr]]/f[rr], {rr, rB, Infinity}, PrecisionGoal -> 12, MaxRecursion -> 30]];
  psiA = ampA; dpsiA = I Qin[rA] ampA; psiB = ampB; dpsiB = I Qup[rB] ampB;
  rsB = rstar[rB];
  eqs = {Psi''[rr] + (w^2 - V[r[rr]]) Psi[rr] == 0, r'[rr] == f[r[rr]]};
  solIn = NDSolve[Join[eqs, {Psi[rsA] == psiA, Psi'[rsA] == dpsiA, r[rsA] == rA}], {Psi, r}, {rr, rsA, rstar[3]},
     PrecisionGoal -> pg, AccuracyGoal -> pg, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}, InterpolationOrder -> All][[1]];
  solUp = NDSolve[Join[eqs, {Psi[rsB] == psiB, Psi'[rsB] == dpsiB, r[rsB] == rB}], {Psi, r}, {rr, rstar[3], rsB},
     PrecisionGoal -> pg, AccuracyGoal -> pg, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}, InterpolationOrder -> All][[1]];
  PsiIn = Psi[rstar[3]] /. solIn; dPsiIn = Psi'[rstar[3]] /. solIn;
  PsiUp = Psi[rstar[3]] /. solUp; dPsiUp = Psi'[rstar[3]] /. solUp;
  W = PsiIn dPsiUp - PsiUp dPsiIn; Ainc = W/(2 I w);
  J = jumpsZRW[l, m, 3, 1, N[bLR, 30], Conjugate[YLRnum[l, m]], Conjugate[dYLRnum[l, m]]];
  P = N[J[[1]]]; dP = N[f[3] J[[2]]];
  ZI = (PsiIn dP - P dPsiIn)/W; ZH = (PsiUp dP - P dPsiUp)/W;
  <|"l" -> l, "m" -> m, "Flux" -> fluxFromZ[l, m, w, ZI, ZH], "u0" -> PsiIn/Ainc, "Du0" -> dPsiIn/Ainc, "v0" -> PsiUp/Ainc, "Dv0" -> dPsiUp/Ainc, "P" -> P, "dP" -> dP|>];

(* small-l test *)
md = photonMode[5, 5]; Print["l = 5, m = 5: Edot_I = ", md["Flux"]["I"], "  Edot_H = ", md["Flux"]["H"]];

(* ::Subsubsection:: *)
(*Toolkit cross-check at l <= 5 (commented out)*)

(* ::Text:: *)
(*This cell needs the Black Hole Perturbation Toolkit (ReggeWheeler and Teukolsky packages) and is commented out. It checks the integrator against the MST solutions of the ReggeWheeler package at l = 5 (the two agree to 10^-7, Sec. III D), and the null source in the Teukolsky equation (last paragraph of Sec. II): with the photon orbit, Upsilon_t = r0^2 p^t = 27, the Teukolsky fluxes coincide with the Zerilli / Regge-Wheeler ones to the precision of the Toolkit's solutions, between nine and twenty significant digits for l <= 5 (the cell runs l = m = 5).*)

(* << ReggeWheeler`
Module[{R, J, PsiIn, dPsiIn, PsiUp, dPsiUp, W, P, dP, ZI, ZH, w = N[5 OmegaLR, 40]},
  R = ReggeWheelerRadial[2, 5, w, "Potential" -> "Zerilli"];
  J = jumpsZRW[5, 5, 3, 1, bLR, Conjugate[YLR[5, 5]], 0]; P = J[[1]]; dP = J[[2]];
  PsiIn = R["In"][3]; dPsiIn = R["In"]'[3]; PsiUp = R["Up"][3]; dPsiUp = R["Up"]'[3];
  W = PsiIn dPsiUp - PsiUp dPsiIn;
  ZI = (PsiIn dP - P dPsiIn)/W; ZH = (PsiUp dP - P dPsiUp)/W;
  Print["Toolkit (MST) l = 5, m = 5: ", fluxFromZ[5, 5, w, ZI, ZH], "   integrator: ", photonMode[5, 5]["Flux"]]];
<< Teukolsky`
photonOrbit = <|"a" -> 0, "p" -> 3, "e" -> 0, "Inclination" -> 1, "Energy" -> 1, "AngularMomentum" -> 3 Sqrt[3],
   "Frequencies" -> <|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)" -> 27|>|>;   (* Upsilon_t = r0^2 p^t = 27 for the photon; the key is the one used by the Toolkit *)
teukPhotonMode[l_, m_, wp_: 40] := Module[{w = N[m OmegaLR, wp], R, S, Z, assoc},
  R = TeukolskyRadial[-2, l, m, 0, w, WorkingPrecision -> wp];
  S = SpinWeightedSpheroidalHarmonicS[-2, l, m, 0];
  Z = Teukolsky`ConvolveSource`Private`ConvolveSourcePointParticleCircular[-2, R, S, photonOrbit];
  assoc = <|"s" -> -2, "l" -> l, "m" -> m, "a" -> 0, "\[Omega]" -> w, "Eigenvalue" -> R["In"]["Eigenvalue"], "Amplitudes" -> Z,
     "Type" -> {"PointParticleCircular"}, "RadialFunctions" -> R, "AngularFunction" -> S, "Method" -> "MST"|>;
  TeukolskyMode[assoc]["Fluxes"]["Energy"]];
Print["Teukolsky null source, l = 5, m = 5: ", teukPhotonMode[5, 5]]; *)

(* ::Subsubsection:: *)
(*Fluxes for l = 10-800 and Fig. 1*)

(* ::Text:: *)
(*Fluxes for a range of l, summed over the five dominant m (j <= 4) and over +-m, and multiplied by l: the first columns of Fig. 1. The paper uses l up to 12800 (stored run, next cells); here the list stops at 800. The mean of the two fluxes sits on kappa already at l of order 100, while the two fluxes approach it as 1 +- sigma_bar/Sqrt[l]. The individual modes are kept in photonModes for the cells that follow.*)

lList = {10, 20, 50, 100, 200, 400, 800};
photonModes = Table[photonMode[l, l - j], {l, lList}, {j, 0, 4}];   (* the five dominant m for each l *)
photonData = Table[{lList[[i]], 2 lList[[i]] Total[#["Flux"]["I"] & /@ photonModes[[i]]], 2 lList[[i]] Total[#["Flux"]["H"] & /@ photonModes[[i]]]}, {i, Length[lList]}];
TableForm[Map[{#[[1]], #[[2]], #[[3]], (#[[2]] + #[[3]])/2, (#[[2]] + #[[3]])/2/kappaTot} &, photonData],
  TableHeadings -> {None, {"l", "l Edot_I", "l Edot_H", "mean", "mean/kappa"}}]

(* Fig. 1 with the data just computed: symbols, the analytic kappa (dashed) and kappa (1 +- sigma_bar/Sqrt[l]) (dotted) *)
Show[
  ListLogLinearPlot[{photonData[[All, {1, 2}]], photonData[[All, {1, 3}]], {#[[1]], (#[[2]] + #[[3]])/2} & /@ photonData},
    PlotMarkers -> {"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]"}, PlotLegends -> {"infinity", "horizon", "mean"}, Joined -> True],
  LogLinearPlot[{kappaTot, kappaTot (1 + sigBar/Sqrt[l]), kappaTot (1 - sigBar/Sqrt[l])}, {l, 8, 1.2 10^4}, PlotStyle -> {Dashed, Dotted, Dotted}],
  Frame -> True, FrameLabel -> {"l", "l Edot_l M^2/E^2"}, PlotRange -> {0.04, 0.1}]

(* ::Subsubsection:: *)
(*The stored run to l = 12800*)

(* ::Text:: *)
(*The stored run up to l = 12800 (photon_big_results.m in the repository root, produced with the same integrator; rows {l, j, Edot_I, Edot_H, u(0), d_r* u(0)}, single m). Reproduces the numbers of Sec. III D: the three-term fit a0 + a2/l + a3 l^(-3/2) of the mean flux (the l^(-1/2) term cancels in the mean) over 400 <= l <= 12800, which gives Sum_j kappa_j = 0.0318263 against the analytic 0.0318265; the individual kappa_j at l = 12800, reproduced to better than 10^-3 in relative terms for j <= 4, with kappa_0 to 10^-4 and kappa_1 to 2 10^-5; and Fig. 1 over the full range (caption of Fig. 1: the five dominant m capture all but 2 10^-6 of the sum, as computed in Section III C).*)

bigRun = Get[FileNameJoin[{repoDir, "photon_big_results.m"}]];
bigL = Union[bigRun[[All, 1]]];
bigData = Table[With[{rows = Select[bigRun, #[[1]] == l &]}, {l, 2 l Total[rows[[All, 3]]], 2 l Total[rows[[All, 4]]]}], {l, bigL}];
fitBig = Fit[Map[{#[[1]], (#[[2]] + #[[3]])/4} &, Select[bigData, #[[1]] >= 400 &]], {1, 1/x, x^(-3/2)}, x];
Print["three-term fit of the single-m mean flux over 400 <= l <= 12800: Sum_j kappa_j = ", fitBig /. x -> Infinity, "   (analytic ", sumKappa, ")"];
ratios12800 = Table[With[{r = First[Select[bigRun, #[[1]] == 12800 && #[[2]] == j &]]}, 12800 (r[[3]] + r[[4]])/2/N[kappaJ[j]]], {j, 0, 4}];
Print["l = 12800, l (Edot_I + Edot_H)/2 / kappa_j for j = 0..4: ", ratios12800];
Print["   relative deviations |ratio - 1|: ", Abs[ratios12800 - 1], ";  largest for j <= 4: ", Max[Abs[ratios12800 - 1]], " (paper: better than 10^-3; kappa_0 to 10^-4, kappa_1 to 2 10^-5)"];
TableForm[Map[{#[[1]], #[[2]], #[[3]], (#[[2]] + #[[3]])/2/kappaTot, Sqrt[#[[1]]] (#[[2]] - #[[3]])/(#[[2]] + #[[3]])} &, bigData], TableHeadings -> {None, {"l", "l Edot_I", "l Edot_H", "mean/kappa", "Sqrt[l](I-H)/(I+H)"}}]
Show[
  ListLogLinearPlot[{bigData[[All, {1, 2}]], bigData[[All, {1, 3}]], {#[[1]], (#[[2]] + #[[3]])/2} & /@ bigData},
    PlotMarkers -> {"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]"}, PlotLegends -> {"infinity", "horizon", "mean"}, Joined -> True],
  LogLinearPlot[{kappaTot, kappaTot (1 + sigBar/Sqrt[l]), kappaTot (1 - sigBar/Sqrt[l])}, {l, 8, 1.5 10^4}, PlotStyle -> {Dashed, Dotted, Dotted}],
  Frame -> True, FrameLabel -> {"l", "l Edot_l M^2/E^2"}, PlotRange -> {0.04, 0.1}]

(* ::Subsubsection:: *)
(*The two contributions to sigma_0 at l = 6400*)

(* ::Text:: *)
(*The two contributions to sigma_0 isolated from the numerical solutions (end of Sec. III D). Since photonMode returns u(0), d_r* u(0), v(0) and d_r* v(0) separately, the asymmetry Sqrt[l] (Edot_I - Edot_H)/(Edot_I + Edot_H) of the m = l mode can be split into the cross term and the barrier asymmetry: at l = 6400, the largest l for which v was computed as well (the stored run only keeps u), the cross term is 1.0389 (analytic 1.0390) and the barrier asymmetry -0.231. The cell also lists the individual kappa_j at the largest l of the list above (l = 800).*)

Print["l = ", Max[lList], ", l (Edot_I + Edot_H)/2 / kappa_j for j = 0..4: ", Table[With[{md = photonModes[[-1, j + 1]]}, Max[lList] (md["Flux"]["I"] + md["Flux"]["H"])/2/N[kappaJ[j]]], {j, 0, 4}]];
Do[md = If[MemberQ[lList, l], photonModes[[Position[lList, l][[1, 1]], 1]], photonMode[l, l]];
  With[{u0 = md["u0"], Du0 = md["Du0"], v0 = md["v0"], Dv0 = md["Dv0"], P = md["P"], dP = md["dP"]},
    Module[{I0 = Abs[u0 dP - P Du0]^2, H0 = Abs[v0 dP - P Dv0]^2, den = (Abs[u0]^2 + Abs[v0]^2) dP^2},
      Print["l = ", l, ":  Sqrt[l](I-H)/(I+H) = ", Sqrt[l] (I0 - H0)/(I0 + H0),
        "   barrier asymmetry = ", Sqrt[l] (Abs[u0]^2 - Abs[v0]^2) dP^2/den,
        "   cross term = ", Sqrt[l] (-2 (Re[u0 Conjugate[Du0]] - Re[v0 Conjugate[Dv0]]) dP P)/den]]],
  {l, {100, 400, 1600, 6400}}];

(* ::Subsubsection:: *)
(*Asymmetry coefficients from the numerical fluxes*)

(* ::Text:: *)
(*The asymmetry coefficients from the numerical fluxes: Sqrt[l] (Edot_I - Edot_H)/(Edot_I + Edot_H) for the individual j (positive for j = 0, negative for j >= 1, against the analytic sigma_j of Section III C) and summed over m (against sigma_bar = 0.61); and the coefficient of the O(1/l) correction to the mean flux, from a fit.*)

TableForm[Table[Prepend[Table[With[{md = photonModes[[Position[lList, l][[1, 1]], j + 1]]}, Sqrt[l] (md["Flux"]["I"] - md["Flux"]["H"])/(md["Flux"]["I"] + md["Flux"]["H"])], {l, {400, Max[lList]}}], sigmaJ[[j + 1, 2]] + sigmaJ[[j + 1, 3]]], {j, 0, 4}],
  TableHeadings -> {Table["j = " <> ToString[j], {j, 0, 4}], {"analytic sigma_j", "l = 400", "l = " <> ToString[Max[lList]]}}]
TableForm[Map[{#[[1]], Sqrt[#[[1]]] (#[[2]] - #[[3]])/(#[[2]] + #[[3]])} &, photonData], TableHeadings -> {None, {"l", "Sqrt[l] (I - H)/(I + H), summed over m"}}]
Print["analytic sigma_bar = ", sigBar, ";  1/l coefficient of mean/kappa - 1 from a fit over l >= 100: ", Fit[Map[{1/#[[1]], (#[[2]] + #[[3]])/2/kappaTot - 1} &, Select[photonData, #[[1]] >= 100 &]], {x}, x] /. x -> 1];

(* ::Section:: *)
(*IV. Extension to Kerr black holes (Table I, Fig. 2)*)

(* ::Subsubsection:: *)
(*Eq. (15): the Kerr light ring*)

(* ::Text:: *)
(*The equatorial light ring of Kerr for signed spin a (a < 0: retrograde orbit; b > 0 and m > 0 throughout). lrData returns the radius r0 = 2 (1 + Cos[2/3 ArcCos[-a]]), the impact parameter b and Upsilon_t from the radial potential R~ = P^2 - Delta (b - a)^2 of equatorial null geodesics, plus the quantities derived from them. The table checks R~(r0) = R~'(r0) = 0 and the two identities R~''(r0) = 6 r0^2 and b^2 - a^2 = 3 r0^2 (beta_b = Sqrt[3] r0) for several spins; the last line proves them symbolically for every light ring, with r0 as the parameter: a = Sqrt[r0] (3 - r0)/2 and b - a = 2 r0 Sqrt[Delta_0]/(r0 - 1), i.e. Eq. (15).*)

s = -2;
Dl[a_, r_] := r^2 - 2 r + a^2;
rpK[a_] := 1 + Sqrt[1 - a^2]; rmK[a_] := 1 - Sqrt[1 - a^2];
rstarK[a_, r_] := r + (2 rpK[a]/(rpK[a] - rmK[a])) Log[(r - rpK[a])/2] - If[a == 0, 0, (2 rmK[a]/(rpK[a] - rmK[a])) Log[(r - rmK[a])/2]];
drdrs[a_, r_] := Dl[a, r]/(r^2 + a^2);
lrData[a_, prec_: 40] := Module[{r0, D0, P0, b, Ups},
  r0 = 2 (1 + Cos[2/3 ArcCos[-a]]); D0 = Dl[a, r0]; P0 = 2 r0 D0/(r0 - 1);
  b = If[a == 0, 3 Sqrt[3], (r0^2 + a^2 - P0)/a];
  Ups = (r0^2 + a^2) P0/D0 + a (b - a);
  N[<|"a" -> a, "r0" -> r0, "b" -> b, "Omega" -> 1/b, "Ups" -> Ups, "D0" -> D0, "P0" -> P0, "betab" -> Sqrt[b^2 - a^2]|>, prec]];
Rtil[a_, b_, r_] := ((r^2 + a^2) - a b)^2 - Dl[a, r] (b - a)^2;
TableForm[Table[With[{d = lrData[a]}, {a, d["r0"], d["b"], d["Ups"], Rtil[a, d["b"], d["r0"]], D[Rtil[a, d["b"], r], r] /. r -> d["r0"],
   D[Rtil[a, d["b"], r], {r, 2}] /. r -> d["r0"], 6 d["r0"]^2, d["b"]^2 - a^2, 3 d["r0"]^2}], {a, {0, 1/2, 9/10, -1/2, -9/10}}] // N,
  TableHeadings -> {None, {"a", "r0", "b", "Ups_t", "R~(r0)", "R~'(r0)", "R~''(r0)", "6 r0^2", "b^2 - a^2", "3 r0^2"}}]
(* symbolic proof: parametrize by r0, with a^2 = r0 (r0-3)^2/4, b - a = 2 r0 Sqrt[Dl]/(r0-1) *)
Clear[r0, a, b];
With[{aa = Sqrt[r0] (3 - r0)/2}, With[{bb = aa + 2 r0 Sqrt[Dl[aa, r0]]/(r0 - 1)},
  Print["for every light ring (1 < r0 < 4): R~''(r0) - 6 r0^2 = ", Simplify[(D[Rtil[aa, bb, r], {r, 2}] /. r -> r0) - 6 r0^2, 1 < r0 < 4], ",  b^2 - a^2 - 3 r0^2 = ", Simplify[bb^2 - aa^2 - 3 r0^2, 1 < r0 < 4]]]]

(* ::Subsubsection:: *)
(*Eq. (17): Omega_theta = lambda_L*)

(* ::Text:: *)
(*Eq. (17): the vertical epicyclic frequency Omega_theta = beta_b/Upsilon_t, the Lyapunov exponent lambda_L = Sqrt[R~''(r0)/(2 Upsilon_t^2)] and the barrier parameter eta_j = (j + 1/2) Omega_theta/lambda_L of the Kerr light rings: Omega_theta = lambda_L, so eta_0 = 1/2 for every spin.*)

kerrIndex[a_] := With[{d = lrData[a]}, With[{Rpp = D[Rtil[a, d["b"], r], {r, 2}] /. r -> d["r0"]},
   {a, N[d["betab"]/d["Ups"], 10], N[Sqrt[Rpp/(2 d["Ups"]^2)], 10], N[(1/2) (d["betab"]/d["Ups"])/Sqrt[Rpp/(2 d["Ups"]^2)], 10]}]];
TableForm[kerrIndex /@ {0, 1/2, 9/10, 99/100, -1/2, -9/10, -99/100}, TableHeadings -> {None, {"a", "Omega_theta", "lambda_L", "eta_0"}}]

(* ::Subsubsection:: *)
(*No superradiance*)

(* ::Text:: *)
(*No superradiance for the resonant modes (remark after Eq. 19): Omega = 1/b exceeds the horizon angular velocity Omega_H = a/(2 r_+) for every Kerr light ring.*)

TableForm[Table[With[{d = lrData[a]}, {N[a], N[d["Omega"], 6], N[a/(2 rpK[a]), 6], N[d["Omega"] - a/(2 rpK[a]), 6]}], {a, {-999/1000, -9/10, -1/2, 0, 1/2, 9/10, 99/100, 999/1000, 9999/10000}}],
  TableHeadings -> {None, {"a", "Omega", "Omega_H", "Omega - Omega_H"}}]

(* ::Subsubsection:: *)
(*Eq. (16): the Teukolsky equation in Schroedinger form, and Im Q*)

(* ::Text:: *)
(*Eq. (16): the Teukolsky radial equation in Schroedinger form, Y'' + Q Y = 0 with R = Delta^(-s/2) (r^2+a^2)^(-1/2) Y (derivatives in r* ), checked by substitution into the Teukolsky equation. Im Q for the resonant modes omega = m/b is s omega [-2 (r-M) P + 4 r Delta]/(r^2+a^2)^2, and the two light-ring conditions make it vanish at r0 for every spin. Near the light ring Im Q = s omega c1 x + O(x^2) is odd in x = r* - r*(r0), i.e. O(omega x/M^2), while the curvature term of Re Q is -(1/2) omega^2 k~ x^2 = O(omega^2 x^2/M^4): in the variable z = (2k)^(1/4) x of Weber's equation, k = omega^2 k~, the equation divided by Sqrt[2k] acquires the odd perturbation I s c1 omega^(-1/2) z/(2k~)^(3/4) = O(l^(-1/2)) z, like the cubic term of the barrier. Far from the light ring the damping q = Im Q/(2p) tends to s/r, the origin of the r^(-2s-1) growth of the Teukolsky function at infinity (checked in the closed-form cell below).*)

Clear[r, a, w, m, lam, s];
Kf = (r^2 + a^2) w - a m;
teuk[R_] := Dl[a, r]^(-s) D[Dl[a, r]^(s + 1) D[R, r], r] + ((Kf^2 - 2 I s (r - 1) Kf)/Dl[a, r] + 4 I s w r - lam) R;
Gf = s (r - 1)/(r^2 + a^2) + r Dl[a, r]/(r^2 + a^2)^2;
DstK[e_] := Dl[a, r]/(r^2 + a^2) D[e, r];
Qexpr = (Kf^2 - 2 I s (r - 1) Kf + Dl[a, r] (4 I s w r - lam))/(r^2 + a^2)^2 - Gf^2 - DstK[Gf];   (* Eq. (16) *)
yppRule = Solve[DstK[DstK[Y[r]]] + Qexpr Y[r] == 0, Y''[r]][[1]];
Print["Teukolsky equation evaluated on R = Dl^(-s/2) (r^2+a^2)^(-1/2) Y, with Y'' + Q Y = 0: ", Simplify[teuk[Dl[a, r]^(-s/2) (r^2 + a^2)^(-1/2) Y[r]] /. yppRule]];
(* Im Q at the light ring for omega = m/b: proportional to -2(r-1)P + 4 r Dl, which vanishes there *)
ImQ = Simplify[ComplexExpand[Im[Qexpr /. m -> w b], TargetFunctions -> {Re, Im}], Assumptions -> {r > 0, a \[Element] Reals, b > 0, w > 0, lam \[Element] Reals}];
Print["Im Q = ", ImQ, "   difference from s w (-2 (r-1) P + 4 r Dl)/(r^2+a^2)^2 at s = -2: ", Simplify[ImQ - s w (-2 (r - 1) ((r^2 + a^2) - a b) + 4 r Dl[a, r])/(r^2 + a^2)^2 /. s -> -2]];
Print["-2 (r0-1) P0 + 4 r0 Dl_0 for a = 0, 1/2, 9/10, -9/10: ", Table[With[{d = lrData[aa]}, -2 (d["r0"] - 1) (d["r0"]^2 + aa^2 - aa d["b"]) + 4 d["r0"] d["D0"]], {aa, {0, 1/2, 9/10, -9/10}}]];
(* the O(l^(-1/2)) scaling: c1 = d/dr* {[-2 (r-1) P + 4 r Dl]/(r^2+a^2)^2} at r0, so that Im Q = s w c1 x + O(x^2); in Schwarzschild Im Q = 2 s w (r - 3)/r^2 exactly *)
Print["Im Q at a = 0: ", Simplify[ImQ /. a -> 0]];
imQc1[aa_] := With[{d = lrData[aa]}, (drdrs[aa, r] D[(-2 (r - 1) ((r^2 + aa^2) - aa d["b"]) + 4 r Dl[aa, r])/(r^2 + aa^2)^2, r]) /. r -> d["r0"]];
ktilde[d_] := Module[{rr, g, a = d["a"]}, g = Rtil[a, d["b"], rr]/(rr^2 + a^2)^2;
  (Dl[a, rr]/(rr^2 + a^2) D[Dl[a, rr]/(rr^2 + a^2) D[g, rr], rr]) /. rr -> d["r0"]];   (* k~ = k/omega^2, the curvature of the barrier *)
Print["c1 = Im Q/(s omega x) at the light ring, a = 0, 1/2, 9/10, -9/10: ", N[imQc1 /@ {0, 1/2, 9/10, -9/10}, 6], "  (a = 0: ", imQc1[0], ")"];
Print["coefficient of the odd perturbation in Weber's equation, s c1/(2 k~)^(3/4), for a = 0, 1/2, 9/10, -9/10: ", N[Table[s imQc1[aa]/(2 ktilde[lrData[aa]])^(3/4), {aa, {0, 1/2, 9/10, -9/10}}] /. s -> -2, 6], " times I omega^(-1/2) z,  with omega^(-1/2) = ", Simplify[Sqrt[bLR/(l - j)], l > j], " in Schwarzschild"];

(* ::Subsubsection:: *)
(*Re Q with the eikonal separation constant*)

(* ::Text:: *)
(*The real part of Q for the resonant modes with the eikonal eigenvalue Lambda = (m - a omega)^2 + (2j+1) omega beta_b + O(1): Re Q = omega^2 R~/(r^2+a^2)^2 - (2j+1) omega beta_b Delta/(r^2+a^2)^2 + O(1), so that the barrier top is at the light ring and its curvature is k = omega^2 k~ (text before Eq. 17). The cell checks the identity.*)

Clear[r, a, w, m, Lam, j, bb, lam0, b];
Print["Re Q - [w^2 R~ - (2j+1) w beta_b Dl]/(r^2+a^2)^2 with the eikonal Lambda, up to the O(1) remainder lam0: ",
  Simplify[((((r^2 + a^2) w - a m)^2 - Dl[a, r] Lam)/(r^2 + a^2)^2 /. {m -> w b, Lam -> (w b - a w)^2 + (2 j + 1) w bb + lam0}) - (w^2 Rtil[a, b, r] - (2 j + 1) w bb Dl[a, r])/(r^2 + a^2)^2]];

(* ::Subsubsection:: *)
(*Toolkit check of the separation constant (commented out)*)

(* ::Text:: *)
(*This cell needs the SpinWeightedSpheroidalHarmonics package of the Toolkit and is commented out. It checks the eigenvalue expansion Lambda = (m - a omega)^2 + (2j+1) omega beta_b + O(1) used above, for m = l - j and omega = m/b (the difference printed must stay O(1) as l grows).*)

(* << SpinWeightedSpheroidalHarmonics`
Table[With[{d = lrData[a], m = l - j}, With[{w = m d["Omega"]},
   {a, l, j, SpinWeightedSpheroidalEigenvalue[-2, l, m, a w] - ((m - a w)^2 + (2 j + 1) w d["betab"])}]],
  {a, {1/2, 9/10}}, {l, {20, 40, 80}}, {j, 0, 2}] // Flatten[#, 2] & // TableForm *)

(* ::Subsubsection:: *)
(*Toolkit check of the spheroidal harmonics at the equator (commented out)*)

(* ::Text:: *)
(*This cell needs the SpinWeightedSpheroidalHarmonics package of the Toolkit and is commented out. It checks the large-omega behaviour of the spheroidal harmonics at the equator that enters the source expansion below: S'(Pi/2)/S(Pi/2) -> (2j+1) s (a+b)/beta_b for even j, S(Pi/2)/S'(Pi/2) -> -s (a+b)/(omega beta_b^2) for odd j, and the normalizations |S(Pi/2)|^2 -> c_j Sqrt[omega beta_b]/(2 Pi^(3/2)) (quoted after Eq. 19; at a = 0, where omega beta_b = m = l - j, it reduces to the spherical result of Sec. III C at leading order in l) and |S'(Pi/2)|^2 -> d_j (omega beta_b)^(3/2)/Pi^(3/2). The ratios printed tend to 1 with O(1/l) corrections.*)

(* << SpinWeightedSpheroidalHarmonics`
swshCheck[a_, l_, j_] := Module[{d = lrData[a], m = l - j, w, g, S, S0, dS0, bb},
  w = m d["Omega"]; g = N[a w]; bb = N[d["betab"]];
  S = Quiet[SpinWeightedSpheroidalHarmonicS[-2, l, m, g, Method -> {"SphericalExpansion", "NumTerms" -> 120}], Eigensystem::arh];
  S0 = S[Pi/2, 0]; dS0 = Derivative[1, 0][S][Pi/2, 0];
  If[EvenQ[j], {a, l, j, (dS0/S0)/((2 j + 1) (-2) (a + d["b"])/bb), Abs[S0]^2/(cj[j] Sqrt[w bb]/(2 Pi^(3/2)))},
               {a, l, j, (S0/dS0)/(2 (a + d["b"])/(w bb^2)), Abs[dS0]^2/(dj[j] (w bb)^(3/2)/Pi^(3/2))}]];   (* s = -2 *)
Table[swshCheck[a, l, j], {a, {1/2, 9/10}}, {l, {50, 100, 200}}, {j, 0, 1}] // Flatten[#, 2] & // N // TableForm *)

(* ::Subsubsection:: *)
(*The source expansion: Ahat_j and Bhat*)

(* ::Text:: *)
(*The source: large-omega expansion of the Teukolsky source projection alpha (the circular-orbit formula of Hughes 2000 in the conventions of the Teukolsky package of the Toolkit) for the null source. The unknowns are the wave at the light ring, Y(0), and its derivative Y'(0) = omega^(1/2) y1; S0 = S(Pi/2), and S'(Pi/2) is related to S0 by the large-omega behaviour checked above (S'' from the angular equation). The O(omega^2) and O(omega^(3/2)) coefficients must vanish (the two cancellations of Sec. IV, which follow from the null condition R~(r0) = 0); the O(omega) coefficient of S0 Y(0) is Ahat_j = (2j+1) Ahat_0 for even j, and the O(omega^(1/2)) coefficient of S0 y1 is Bhat for odd j. In Schwarzschild Ahat_0 = 9 Sqrt[27]/4 = 3 r0 b/4 and Bhat = -I Ahat_0. The table lists the two leading coefficients, Ahat_j/Ahat_0 and Bhat for seven spins and j <= 3.*)

alphaExpand[a_, j_] := Module[{d = lrData[a, 40], r0, b, bb, D0, P0, w, u, SS, YY0, yy1, lam0, m, c, lam, A, Kt, L, S0, dS0, d2S0, L2S, L1L2S, rho, rhob, Sig,
    Ann0, Anmb0, Anmb1, Ambmb0, Ambmb1, Ambmb2, rc, tc, Cnn, Cnmb, Cmbmb, r, pref, pref0, dpref0, R, dR, d2R, alpha, poly, co, sqw},
  {r0, b, bb, D0, P0} = d /@ {"r0", "b", "betab", "D0", "P0"};
  m = w b; c = a w; lam = (m - c)^2 + (2 j + 1) w bb + lam0; A = lam - c^2 + 2 m c; Kt = w P0; L = -m + c;
  If[EvenQ[j], S0 = SS; dS0 = (2 j + 1) s (a + b)/bb SS, dS0 = SS; S0 = -s (a + b)/bb^2 SS/w];
  d2S0 = (m^2 - s - A) S0;
  L2S = dS0 + L S0; L1L2S = d2S0 + 2 L dS0 - 2 S0 + L^2 S0;
  rho = -1/r0; rhob = -1/r0; Sig = r0^2;
  Ann0 = -rho^-2 rhob^-1 (Sqrt[2] D0)^-2 (rho^-1 L1L2S + 2 I a L2S);
  Anmb0 = rho^-3 (Sqrt[2] D0)^-1 ((2 rho - I Kt/D0) L2S);
  Anmb1 = -rho^-3 (Sqrt[2] D0)^-1 L2S;
  Ambmb0 = (Kt^2 S0 rhob)/(4 D0^2 rho^3) + (I Kt S0 (1 - r0 + D0 rho) rhob)/(2 D0^2 rho^3) + (I r0 S0 rhob w)/(2 D0 rho^3);
  Ambmb1 = -rho^-3 rhob S0/2 (I Kt/D0 - rho);
  Ambmb2 = -rho^-3 rhob S0/4;
  rc = P0/(2 Sig); tc = I (b - a)/(Sqrt[2] r0);
  {Cnn, Cnmb, Cmbmb} = {rc^2, rc tc, tc^2};
  pref = (r^2 - 2 r + a^2) (r^2 + a^2)^(-1/2);
  pref0 = pref /. r -> r0; dpref0 = D[pref, r] /. r -> r0;
  R = pref0 YY0; dR = dpref0 YY0 + pref0 ((r0^2 + a^2)/D0) sqw yy1;   (* sqw stands for Sqrt[w] *)
  d2R = (-(-lam + 2 I r0 s 2 w + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) w) + (-a m + (a^2 + r0^2) w)^2)/D0) R - (-2 + 2 r0) (1 + s) dR)/D0;
  alpha = (Ann0 Cnn + Anmb0 Cnmb + Ambmb0 Cmbmb) R - (Anmb1 Cnmb + Ambmb1 Cmbmb) dR + Ambmb2 Cmbmb d2R;
  poly = Expand[(alpha /. {w -> u^2, sqw -> u}) u^8];   (* Laurent polynomial in u = Sqrt[omega], shifted by u^8 *)
  co[p_] := Chop[Coefficient[poly, u, 2 p + 8], 10^-25];
  <|"data" -> d, "omega2" -> co[2], "omega32" -> co[3/2], "Ahat" -> Coefficient[co[1], SS YY0], "Bhat" -> Coefficient[co[1/2], SS yy1], "omega1" -> co[1], "omega12" -> co[1/2]|>];
s = -2;
ae = alphaExpand[0, 0];
Print["a = 0, j = 0: omega^2 coefficient = ", ae["omega2"], ",  omega^(3/2) coefficient = ", ae["omega32"]];
Print["   Ahat_0 = ", ae["Ahat"], "  (expected 9 Sqrt[27]/4 = ", N[9 Sqrt[27]/4, 20], "),  Bhat = ", alphaExpand[0, 1]["Bhat"], "  (expected -I Ahat_0)"];
(* columns: a, j, omega^2 and omega^(3/2) coefficients (must vanish), Ahat_j/Ahat_0 (= 2j+1 for even j, 0 for odd j), and for odd j the coefficient Bhat
   (the O(omega^(1/2)) coefficient returned for even j is a different, subleading quantity and is not used) *)
TableForm[Table[With[{e = alphaExpand[a, j]}, {a, j, Chop[e["omega2"]], Chop[e["omega32"]], e["Ahat"]/alphaExpand[a, 0]["Ahat"], If[OddQ[j], e["Bhat"], "-"]}], {a, {0, 1/2, 9/10, 99/100, -1/2, -9/10, -99/100}}, {j, 0, 3}] // Flatten[#, 1] &,
  TableHeadings -> {None, {"a", "j", "omega^2 coeff.", "omega^(3/2) coeff.", "Ahat_j/Ahat_0", "Bhat"}}]

(* ::Subsubsection:: *)
(*Eqs. (18)-(20): amplitude integrals, kappa_j(a) and the horizon/infinity ratio*)

(* ::Text:: *)
(*Eqs. (18)-(20) evaluated numerically. The amplitude integrals J_inf and h_inf of Eq. (18) and the factor N = r0^(2s) E^(-2 s h_inf) E^(-2 J_inf), which tracks the WKB amplitude of the incident wave from infinity to the light ring, |Y_inc|^2 = |B_inc|^2 (omega/p) N; the factor omega/p is the usual WKB amplitude factor |Q|^(-1/2) referred to infinity, where p -> omega, and it is the same factor 1/p_B that normalizes the incident wave in the model equation of Section III B. With it, |u(0)|^2 = N (omega/(2k)^(1/4)) (Pi/Sqrt[2]) F(eta_j), and the cell evaluates the barrier curvature k~, kappa_j(a) of Eq. (19) for even j and its odd counterpart, the horizon-side factors (J_H, the Wronskian normalization W = 2 I C_out (r_+^2 + a^2)(varpi - 2 I s epsilon_+) and the Teukolsky-Press factor), whose ratio to the infinity-side ones is 1 for every spin to the working precision (footnote of Sec. IV), and the identity |Bhat/Ahat_0|^2 Sqrt[2 k~] beta_b = 2, which makes kappa_j(a)/kappa_0(a) spin independent. At a = 0 the Teukolsky bookkeeping must reproduce the Zerilli / Regge-Wheeler kappa_j of Section III C, which it does to all digits printed.*)

qdrs[a_, b_, rr_] := s (-2 (rr - 1) ((rr^2 + a^2) - a b) + 4 rr Dl[a, rr])/(2 Dl[a, rr] Sqrt[Rtil[a, b, rr]]);   (* q dr*/dr of Eq. (18) *)
drsdr[a_, rr_] := (rr^2 + a^2)/Dl[a, rr];
Jinf[d_] := NIntegrate[qdrs[d["a"], d["b"], rr] - s/rr drsdr[d["a"], rr], {rr, d["r0"], d["r0"] + 1, d["r0"] + 10, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 15, MaxRecursion -> 20];
hinf[d_] := NIntegrate[(rr^2 + d["a"]^2)/(rr Dl[d["a"], rr]) - 1/rr, {rr, d["r0"], d["r0"] + 10, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 15, MaxRecursion -> 20];
Jhor[d_] := NIntegrate[qdrs[d["a"], d["b"], rr] + s (rr - 1)/(rr^2 + d["a"]^2) drsdr[d["a"], rr], {rr, rpK[d["a"]], d["r0"]}, WorkingPrecision -> 30, PrecisionGoal -> 15, MaxRecursion -> 20];
kappaKerr[a_, j_] := Module[{ae = alphaExpand[a, j], d, r0, b, bb, Ups, D0, kt, J0, JH, hi, N0, rp, OmH, kap, NH, alphaH, x, kI, Du},
  d = ae["data"]; {r0, b, bb, Ups, D0} = d /@ {"r0", "b", "betab", "Ups", "D0"};
  kt = ktilde[d]; J0 = Jinf[d]; JH = Jhor[d]; hi = hinf[d];
  N0 = r0^(2 s) Exp[-2 s hi] Exp[-2 J0];                                   (* N of Eq. (18) *)
  rp = rpK[a]; OmH = a/(2 rp); kap = 1 - b OmH;                           (* kap = varpi/omega *)
  NH = (2 rp) D0^s Exp[-2 JH]/((2 rp)^2 kap^2);                            (* horizon-side amplitude factor and Wronskian normalization *)
  alphaH = 256 (2 rp)^5 kap^5/(b - a)^8;                                  (* Teukolsky-Press factor at leading order *)
  x = j + 1/2;
  kI = If[EvenQ[j],
    Sqrt[2 Pi] b cj[j] Sqrt[bb] N0 Fpar[x] Abs[ae["Ahat"]]^2/((2 kt)^(1/4) Ups^2),                     (* Eq. (19) *)
    Du = (2 kt)^(1/4) N0 Sqrt[2] Pi Exp[-Pi x/2]/(Cosh[Pi x] G14[x]);
    4 Pi Abs[ae["Bhat"]]^2 dj[j] bb^(3/2)/Pi^(3/2) Du b/Ups^2];                                      (* odd counterpart *)
  <|"a" -> a, "j" -> j, "kappa" -> kI, "ratioHI" -> alphaH kap NH/N0, "ktilde" -> kt, "Ahat" -> ae["Ahat"], "N" -> N0, "Ups" -> Ups, "r0" -> r0, "b" -> b, "Jinf" -> J0, "hinf" -> hi, "JH" -> JH|>];
Print["a = 0: J_inf = ", kappaKerr[0, 0]["Jinf"], ",  h_inf = ", kappaKerr[0, 0]["hinf"], ",  N = ", kappaKerr[0, 0]["N"], "  (16/r0^6 = ", N[16/3^6, 10], ")"];
Print["horizon/infinity ratio - 1 for a = 0, 1/2, 9/10, 99/100, -1/2, -9/10, -99/100: ", Table[kappaKerr[a, 0]["ratioHI"] - 1, {a, {0, 1/2, 9/10, 99/100, -1/2, -9/10, -99/100}}], "  (the deviations are the precision of the numerical integrals)"];
(* in Schwarzschild Bhat = -I Ahat_0; for general spin |Bhat/Ahat_0|^2 Sqrt[2 k~] beta_b = 2 *)
TableForm[Table[With[{e0 = alphaExpand[a, 0], e1 = alphaExpand[a, 1], d = lrData[a]}, {a, e1["Bhat"]/e0["Ahat"], Abs[e1["Bhat"]/e0["Ahat"]]^2 Sqrt[2 ktilde[d]] d["betab"]}], {a, {0, 1/2, 9/10, 99/100, -1/2, -9/10, -99/100}}],
  TableHeadings -> {None, {"a", "Bhat/Ahat_0", "|Bhat/Ahat_0|^2 Sqrt[2 k~] beta_b"}}]
(* a = 0 must reproduce the Zerilli / Regge-Wheeler kappa_j of Section III C *)
TableForm[Table[{j, kappaKerr[0, j]["kappa"], N[kappaJ[j], 12], kappaKerr[0, j]["ratioHI"]}, {j, 0, 4}], TableHeadings -> {None, {"j", "kappa_j from Eq. (19), a = 0", "kappa_j of Section III C", "horizon/infinity"}}]

(* ::Subsubsection:: *)
(*Table I*)

(* ::Text:: *)
(*Table I: for each spin, r0, b, Upsilon_t, k~, |Ahat_0|, N, kappa_0(a) and g(a) = kappa_0(a)/kappa_0(0) (Eq. 20), as evaluated numerically; the last column here is the analytic ratio of the horizon to the infinity flux, which is 1 (the last column of the paper's Table I compares with the numerical runs and is reproduced in a later cell). Then the spin independence of kappa_1/kappa_0, the values g = 0.99 at a = 0.98 and 0.62 at a = 0.995 quoted after Eq. (21), and the position of the maximum of g for prograde orbits.*)

spins = {0, 1/2, 9/10, 99/100, -1/2, -9/10, -99/100};
kappa00 = kappaKerr[0, 0]["kappa"];
tableI = Table[With[{k = kappaKerr[a, 0]}, {N[a], N[k["r0"], 5], N[k["b"], 5], N[k["Ups"], 5], N[k["ktilde"], 4], N[Abs[k["Ahat"]], 4], N[k["N"], 4], N[k["kappa"], 5], N[k["kappa"]/kappa00, 4], N[k["ratioHI"], 10]}], {a, spins}];
TableForm[tableI, TableHeadings -> {None, {"a", "r0", "b", "Ups_t", "k~", "|Ahat_0|", "N", "kappa_0(a)", "g(a)", "H/I"}}]
Print["kappa_1(a)/kappa_0(a) for a = 0, 1/2, 9/10: ", Table[kappaKerr[a, 1]["kappa"]/kappaKerr[a, 0]["kappa"], {a, {0, 1/2, 9/10}}]];
Print["g(a) for a = 0.7, 0.8, 0.95, 0.98, 0.995: ", Table[{N[a], N[kappaKerr[a, 0]["kappa"]/kappa00, 5]}, {a, {7/10, 8/10, 95/100, 98/100, 995/1000}}]];
(* location of the maximum of g(a) for prograde orbits *)
gfun[x_?NumericQ] := Re[kappaKerr[Rationalize[x, 0], 0]["kappa"]]/kappa00;
Print["maximum of g(a): ", FindMaximum[gfun[x], {x, 0.83, 0.80, 0.86}], "   (Eq. 21: 6 Sqrt[3] - 9 = ", N[6 Sqrt[3] - 9, 6], " at a = ", N[3^(1/4) (3 - Sqrt[3])/2, 6], ")"];

(* ::Subsubsection:: *)
(*Closed forms and Eq. (21)*)

(* ::Text:: *)
(*Closed forms: Eq. (15), the expressions quoted in the text of Sec. IV (Upsilon_t, beta_b, k~, Ahat_0, Bhat), Eq. (18) with N = 16 M^2/r0^6, and Eq. (21). With the light ring parametrized by its radius r0, a = Sqrt[r0] (3 - r0)/2, every ingredient of Eq. (19) is elementary: R~ = (r - r0)^2 r (r + 2 r0) has a double root at r0, so the amplitude integrals are logarithms. This cell derives b, Upsilon_t, beta_b, k~, J_inf, h_inf, N, Ahat_0, Bhat and g(a) symbolically, checks the two source cancellations and the identity |Bhat/Ahat_0|^2 Sqrt[2 k~] beta_b = 2 for every r0, gives the maximum of g and its value at r0 = 4, and proves that the horizon/infinity ratio is exactly 1: J_H is made rational by the Euler substitution Sqrt[r (r + 2 r0)] = t - r, giving the closed form E^(2 J_H) = [r0^2 + r0 + 4 Sqrt[4 - r0] - 8]^4/r0^8 quoted in the text, and the radicals of r_+ are denested with r0 = 4 - u^2 (u between 1 and Sqrt[3] for prograde light rings, between 0 and 1 for retrograde ones). All simplifications assume 1 < r0 < 4, i.e. both branches.*)

Clear[r0, r, t, u];
aS = Sqrt[r0] (3 - r0)/2; D0S = (r0 - 1)^2 r0/4; P0S = 2 r0 D0S/(r0 - 1);
bS = Simplify[aS + 2 r0 Sqrt[D0S]/(r0 - 1), 1 < r0 < 4];
UpsS = Simplify[(r0^2 + aS^2) P0S/D0S + aS (bS - aS), 1 < r0 < 4];
betabS = Simplify[Sqrt[bS^2 - aS^2], 1 < r0 < 4];
quadS = Factor[Simplify[Rtil[aS, bS, r]/(r - r0)^2, 1 < r0 < 4]];
ktS = Simplify[(Dl[aS, r]/(r^2 + aS^2) D[Dl[aS, r]/(r^2 + aS^2) D[Rtil[aS, bS, r]/(r^2 + aS^2)^2, r], r]) /. r -> r0, 1 < r0 < 4];
Print["Eq. (15): b = ", bS, ",  R~/(r - r0)^2 = ", quadS, ";  Ups_t = ", UpsS, ",  beta_b = ", betabS, ",  k~ = ", ktS, ",  k~ - 3 r0^4 (r0-1)^4/(8 (r0^2+a^2)^4) = ", Simplify[ktS - 3 r0^4 (r0 - 1)^4/(8 (r0^2 + aS^2)^4)]];
sqS = (r - r0) Sqrt[r (r + 2 r0)];                                 (* Sqrt[R~] for r > r0 *)
JinfS = Integrate[Simplify[s (-2 (r - 1) ((r^2 + aS^2) - aS bS) + 4 r Dl[aS, r])/(2 Dl[aS, r] sqS) - s/r (r^2 + aS^2)/Dl[aS, r], 1 < r0 < 4], {r, r0, Infinity}, Assumptions -> 1 < r0 < 4];
hinfS = Integrate[(r^2 + aS^2)/(r Dl[aS, r]) - 1/r, {r, r0, Infinity}, Assumptions -> 1 < r0 < 4];
NS = Simplify[r0^(2 s) Exp[-2 s hinfS] Exp[-2 JinfS], 1 < r0 < 4];
Print["Eq. (18): J_inf = ", Simplify[JinfS, 1 < r0 < 4], ",  h_inf = ", Simplify[hinfS, 1 < r0 < 4], ",  N = ", NS, ",  N - 16/r0^6 = ", Simplify[NS - 16/r0^6]];
Print["r q at infinity (q = Im Q/(2p) -> s/r, text after Eq. 18): ", Limit[r s (-2 (r - 1) ((r^2 + aS^2) - aS bS) + 4 r Dl[aS, r])/(2 Dl[aS, r] sqS) Dl[aS, r]/(r^2 + aS^2), r -> Infinity, Assumptions -> 1 < r0 < 4], ";  omega beta_b at a = 0 (r0 = 3) for omega = m/b: ", Simplify[(m/bS) betabS /. r0 -> 3]];
(* source coefficients: alphaExpand evaluated on the symbolic light-ring data *)
lrS = <|"a" -> aS, "r0" -> r0, "b" -> bS, "Omega" -> 1/bS, "Ups" -> UpsS, "D0" -> D0S, "P0" -> P0S, "betab" -> betabS|>;
{ae0, ae1} = Block[{lrData}, lrData[_, ___] := lrS; {alphaExpand[aS, 0], alphaExpand[aS, 1]}];
Print["omega^2 and omega^(3/2) coefficients of the source, j = 0 and 1: ", Simplify[{ae0["omega2"], ae0["omega32"], ae1["omega2"], ae1["omega32"]}, 1 < r0 < 4]];
AhS = Simplify[ae0["Ahat"], 1 < r0 < 4]; BhS = Simplify[ae1["Bhat"], 1 < r0 < 4];
Print["Ahat_0 = ", AhS, " (difference from Sqrt[3] r0^4/(4 Sqrt[r0^2+a^2]): ", Simplify[AhS - Sqrt[3] r0^4/(4 Sqrt[r0^2 + aS^2]), 1 < r0 < 4], "),  Bhat = ", BhS, " (difference from -I r0^(5/2) Sqrt[r0^2+a^2]/(2 (r0-1)): ", Simplify[BhS + I r0^(5/2) Sqrt[r0^2 + aS^2]/(2 (r0 - 1)), 1 < r0 < 4], ")"];
Print["|Bhat/Ahat_0|^2 Sqrt[2 k~] beta_b = ", Simplify[ComplexExpand[Abs[BhS]^2/Abs[AhS]^2] Sqrt[2 ktS] betabS, 1 < r0 < 4]];
kap0S = bS Sqrt[betabS] NS ComplexExpand[Abs[AhS]^2]/((2 ktS)^(1/4) UpsS^2);        (* kappa_0(a) up to spin-independent factors *)
gS = Simplify[kap0S/(kap0S /. r0 -> 3), 1 < r0 < 4];
Print["Eq. (21): g = ", gS, ";  maximum: ", Simplify[Maximize[{gS, 1 < r0 < 4}, r0]], ",  a/M at the maximum: ", Simplify[aS /. r0 -> Sqrt[3]], " = ", N[aS /. r0 -> Sqrt[3], 6], ",  g(r0 = 4) = ", gS /. r0 -> 4, " = ", N[gS /. r0 -> 4, 4]];
(* horizon side *)
rpS = 1 + Sqrt[1 - aS^2];
integH = Simplify[s (-2 (r - 1) ((r^2 + aS^2) - aS bS) + 4 r Dl[aS, r])/(2 Dl[aS, r] (r0 - r) Sqrt[r (r + 2 r0)]) + s (r - 1)/Dl[aS, r], 1 < r0 < 4];
rT = t^2/(2 (t + r0));
integT = Simplify[((integH /. Sqrt[r (r + 2 r0)] -> t - r) /. r -> rT) D[rT, t], {t > 0, 1 < r0 < 4}];
antiT = Integrate[integT, t];
JHS = (antiT /. t -> r0 + Sqrt[r0 (r0 + 2 r0)]) - (antiT /. t -> rpS + Sqrt[rpS (rpS + 2 r0)]);
kapS = 1 - bS aS/(2 rpS);
ratioHI = 256 (2 rpS)^5 kapS^5/(bS - aS)^8 kapS (2 rpS) D0S^s Exp[-2 JHS]/((2 rpS)^2 kapS^2)/NS;
Print["Exp[2 J_H] r0^8/(r0^2 + r0 + 4 Sqrt[4 - r0] - 8)^4 - 1, denested with r0 = 4 - u^2: prograde ", FullSimplify[(Exp[2 JHS] r0^8/(r0^2 + r0 + 4 Sqrt[4 - r0] - 8)^4 - 1) /. r0 -> 4 - u^2, 1 < u < Sqrt[3]], ", retrograde ", FullSimplify[(Exp[2 JHS] r0^8/(r0^2 + r0 + 4 Sqrt[4 - r0] - 8)^4 - 1) /. r0 -> 4 - u^2, 0 < u < 1]];
Print["horizon/infinity ratio - 1: prograde ", FullSimplify[(ratioHI - 1) /. r0 -> 4 - u^2, 1 < u < Sqrt[3]], ", retrograde ", FullSimplify[(ratioHI - 1) /. r0 -> 4 - u^2, 0 < u < 1], ";  numerically at r0 = 23/10: ", N[ratioHI /. r0 -> 23/10, 30]];

(* ::Subsubsection:: *)
(*g(a) on a dense grid*)

(* ::Text:: *)
(*The spin factor g(a) on a dense grid, from the numerical evaluation of Eq. (19) (solid; the lower panel of Fig. 2) against the closed form of Eq. (21) (dashed).*)

gGrid = Join[Range[-99/100, -90/100, 3/100], Range[-9/10, 9/10, 1/20], Range[9/10, 995/1000, 1/200]];
gCurve = Table[{N[a], N[kappaKerr[a, 0]["kappa"]/kappa00]}, {a, gGrid}];
gClosed[a_] := With[{r = 2 (1 + Cos[2/3 ArcCos[-a]])}, 27 (r - 1)/(r^2 (r + 3))];   (* Eq. (21) *)
Print["max |g_numerical/g_closed - 1| on the grid: ", Max[Abs[gCurve[[All, 2]]/(gClosed /@ N[gGrid]) - 1]]];
Show[ListLinePlot[gCurve, Frame -> True, FrameLabel -> {"a/M", "g(a)"}, PlotRange -> {{-1, 1}, {0.5, 1.5}}, GridLines -> {None, {1}}], Plot[gClosed[a], {a, -1, 1}, PlotStyle -> {Red, Dashed}]]

(* ::Subsubsection:: *)
(*Kerr Teukolsky integrator (commented out)*)

(* ::Text:: *)
(*This cell needs the SpinWeightedSpheroidalHarmonics package of the Toolkit and is commented out. Numerical Teukolsky fluxes for the photon on the Kerr light ring, the data of Fig. 2: direct integration of Y'' + Q Y = 0 with WKB boundary data, spheroidal harmonics from the Toolkit. The 'in' solution is integrated in the Schroedinger form from the horizon; the 'up' solution is obtained from the Sasaki-Nakamura equation, whose potential is short ranged. Conventions of the Teukolsky package of the Toolkit (R_in with unit transmission at the horizon, R_up with unit transmission at infinity). At a = 0 the result must agree with photonMode. The l = 800 data of the paper were produced with the Python implementation of the same algorithm in the kerr_py directory of the repository (eighth-order Dormand-Prince integrator, with its own spectral spheroidal harmonics); the cell after the next one reads those stored results.*)

(* << SpinWeightedSpheroidalHarmonics`
QfunK[a_, m_, w_, lam_] := Module[{r, K, G, Q}, K = (r^2 + a^2) w - a m; G = s (r - 1)/(r^2 + a^2) + r Dl[a, r]/(r^2 + a^2)^2;
  Q = (K^2 - 2 I s (r - 1) K + Dl[a, r] (4 I s w r - lam))/(r^2 + a^2)^2 - G^2 - drdrs[a, r] D[G, r]; Function @@ {r, Q}];
QwkbK[a_, m_, w_, lam_, sign_, niter_: 3] := Module[{r, Q, q}, Q = QfunK[a, m, w, lam][r]; q = sign Sqrt[Q];
  Do[q = sign Sqrt[Q + I drdrs[a, r] D[q, r]], {niter}]; Function @@ {r, q}];
alphaSource[a_, r0_, EE_, LL_, Ut_, m_, w_, lam_, S0_, dS0_, d2S0_, R_, dR_] :=
 Module[{th0 = Pi/2, D0, Kt, d2R, L1, L2, L2S, L2p, L1Sp, L1L2S, rho, rhob, Sig, Ann0, Anmb0, Anmb1, Ambmb0, Ambmb1, Ambmb2, rc, tc, Cnn, Cnmb, Cmbmb},
  D0 = r0^2 - 2 r0 + a^2; Kt = (r0^2 + a^2) w - m a;
  d2R = (-(-lam + 2 I r0 s 2 w + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) w) + (-a m + (a^2 + r0^2) w)^2)/D0) R - (-2 + 2 r0) (1 + s) dR)/D0;
  L1 = -m/Sin[th0] + a w Sin[th0] + Cos[th0]/Sin[th0]; L2 = -m/Sin[th0] + a w Sin[th0] + 2 Cos[th0]/Sin[th0];
  L2S = dS0 + L2 S0; L2p = m Cos[th0]/Sin[th0]^2 + a w Cos[th0] - 2/Sin[th0]^2; L1Sp = d2S0 + L1 dS0; L1L2S = L1Sp + L2p S0 + L2 dS0 + L1 L2 S0;
  rho = -1/(r0 - I a Cos[th0]); rhob = -1/(r0 + I a Cos[th0]); Sig = 1/(rho rhob);
  Ann0 = -rho^-2 rhob^-1 (Sqrt[2] D0)^-2 (rho^-1 L1L2S + 3 I a Sin[th0] L1 S0 + 3 I a Cos[th0] S0 + 2 I a Sin[th0] dS0 - I a Sin[th0] L2 S0);
  Anmb0 = rho^-3 (Sqrt[2] D0)^-1 ((rho + rhob - I Kt/D0) L2S + (rho - rhob) a Sin[th0] Kt/D0 S0);
  Anmb1 = -rho^-3 (Sqrt[2] D0)^-1 (L2S + I (rho - rhob) a Sin[th0] S0);
  Ambmb0 = (Kt^2 S0 rhob)/(4 D0^2 rho^3) + (I Kt S0 (1 - r0 + D0 rho) rhob)/(2 D0^2 rho^3) + (I r0 S0 rhob w)/(2 D0 rho^3);
  Ambmb1 = -rho^-3 rhob S0/2 (I Kt/D0 - rho); Ambmb2 = -rho^-3 rhob S0/4;
  rc = (EE (r0^2 + a^2) - a LL)/(2 Sig); tc = rho (I Sin[th0] (a EE - LL/Sin[th0]^2))/Sqrt[2];
  {Cnn, Cnmb, Cmbmb} = {rc^2, rc tc, tc^2};
  (Ann0 Cnn + Anmb0 Cnmb + Ambmb0 Cmbmb) R - (Anmb1 Cnmb + Ambmb1 Cmbmb) dR + Ambmb2 Cmbmb d2R];
fluxKerr[a_, m_, w_, lam_, ZI_, ZH_] := Module[{rh = rpK[a], OmH, kap, eps, C2, al}, OmH = a/(2 rh); kap = w - m OmH; eps = Sqrt[1 - a^2]/(4 rh);
  C2 = ((lam + 2)^2 + 4 a m w - 4 a^2 w^2) (lam^2 + 36 m a w - 36 a^2 w^2) + (2 lam + 3) (96 a^2 w^2 - 48 m a w) + 144 w^2 (1 - a^2);
  al = 256 (2 rh)^5 kap (kap^2 + 4 eps^2) (kap^2 + 16 eps^2) w^3/C2;
  <|"I" -> Abs[ZI]^2/(4 Pi w^2), "H" -> al Abs[ZH]^2/(4 Pi w^2)|>];
swshData[a_, l_, m_, w_, nterms_: 120] := Module[{g = N[a w], lam, S},
  (* machine-precision spherical expansion; eigenvalue and harmonic at theta = Pi/2 *)
  lam = Quiet[SpinWeightedSpheroidalEigenvalue[-2, l, m, g, Method -> {"SphericalExpansion", "NumTerms" -> nterms}], {Eigensystem::arhm, Eigensystem::arh}];
  If[a == 0,
    N[{lam, SpinWeightedSphericalHarmonicY[-2, l, m, Pi/2, 0], Derivative[0, 0, 0, 1, 0][SpinWeightedSphericalHarmonicY][-2, l, m, Pi/2, 0], Derivative[0, 0, 0, 2, 0][SpinWeightedSphericalHarmonicY][-2, l, m, Pi/2, 0]}],
    S = Quiet[SpinWeightedSpheroidalHarmonicS[-2, l, m, g, Method -> {"SphericalExpansion", "NumTerms" -> nterms}], {Eigensystem::arhm, Eigensystem::arh}];
    N[{lam, S[Pi/2, 0], Derivative[1, 0][S][Pi/2, 0], Derivative[2, 0][S][Pi/2, 0]}]]];
(* "Up" solution from the Sasaki-Nakamura equation (short-ranged potential), Toolkit conventions: R_up -> r^(-2s-1) e^(i w rstar) at infinity.
   X satisfies X_** - F X_* - U X = 0 (derivatives in rstar); Z = X exp(-Int F/2 drstar) satisfies Z'' + Qsn Z = 0, with WKB data at rB. *)
snUpAt[a_, w_, m_, lam_, r0_, rB_: 40] := Module[{r, Dl, f, intQ, intF, Rmax, c0, c1, c2, c3, c4, eta, K, V, beta, alpha, U1, G, F, U, Qsn, q, Ff, Uf, qf, ff, Iq, IF, ZB, dZB, XB, dXB, X, rv, x, rr, sol, X0, dXs0, dX0, d2X0, f0, fp, scale, Xf, chi, Rup, dRup, rules, rsB, rs0},
  Dl = r^2 - 2 r + a^2; f = Dl/(r^2 + a^2);
  c0 = -12 I w + lam (lam + 2) - 12 a w (a w - m);
  c1 = 8 I a (3 a w - lam (a w - m));
  c2 = -24 I a (a w - m) + 12 a^2 (1 - 2 (a w - m)^2);
  c3 = 24 I a^3 (a w - m) - 24 a^2;
  c4 = 12 a^4;
  eta = c0 + c1/r + c2/r^2 + c3/r^3 + c4/r^4;
  K = (r^2 + a^2) w - m a;
  V = -((K^2 + 4 I (r - 1) K)/Dl) + 8 I w r + lam;
  beta = 2 Dl*(-I K + r - 1 - 2 Dl/r);
  alpha = -I K beta/Dl^2 + 3 I D[K, r] + lam + 6 Dl/r^2;
  U1 = V + Dl^2/beta (D[2 alpha + D[beta, r]/Dl, r] - D[eta, r]/eta (alpha + D[beta, r]/Dl));
  G = -((2 (r - 1))/(r^2 + a^2)) + (r Dl)/(r^2 + a^2)^2;
  F = D[eta, r]/eta Dl/(r^2 + a^2);
  U = (Dl U1)/(r^2 + a^2)^2 + G^2 + (Dl D[G, r])/(r^2 + a^2) - F G;
  Qsn = -U - F^2/4 - f D[F, r]/2;
  q = Sqrt[Qsn]; Do[q = Sqrt[Qsn + I f D[q, r]], {2}];
  Ff = Function @@ {r, F}; Uf = Function @@ {r, U}; qf = Function @@ {r, q}; ff = Function @@ {r, f};
  rsB = N[rstarK[a, rB]]; rs0 = N[rstarK[a, r0]];
  (* the WKB expressions are large: evaluate them only at numerical points *)
  intQ[rr_?NumericQ] := I (qf[rr] - w)/ff[rr]; intF[rr_?NumericQ] := (Ff[rr]/2)/ff[rr];
  (* both integrands fall off as 1/r^2: integrate numerically up to Rmax and add the analytic tail c/Rmax, c = integrand(Rmax) Rmax^2 *)
  Rmax = 10^4 rB;
  Iq = NIntegrate[intQ[rr], {rr, rB, 10 rB, 100 rB, 1000 rB, Rmax}, PrecisionGoal -> 9, MaxRecursion -> 30, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}] + intQ[Rmax] Rmax;
  IF = If[a == 0, 0, NIntegrate[intF[rr], {rr, rB, 10 rB, 100 rB, 1000 rB, Rmax}, PrecisionGoal -> 9, MaxRecursion -> 30, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}] + intF[Rmax] Rmax];   (* F = 0 in Schwarzschild *)
  ZB = Exp[I w rsB] Exp[-Iq]; dZB = I qf[rB] ZB; XB = Exp[-IF] ZB; dXB = Exp[-IF] (dZB + Ff[rB]/2 ZB);
  sol = NDSolve[{X''[x] == Ff[rv[x]] X'[x] + Uf[rv[x]] X[x], rv'[x] == ff[rv[x]], X[rsB] == XB, X'[rsB] == dXB, rv[rsB] == N[rB]}, {X, rv}, {x, rs0, rsB},
     PrecisionGoal -> 12, AccuracyGoal -> 12, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9, "StiffnessTest" -> False}, InterpolationOrder -> All][[1]];
  X0 = X[rs0] /. sol; dXs0 = X'[rs0] /. sol;
  f0 = ff[r0]; fp = D[f, r] /. r -> r0;
  dX0 = dXs0/f0; d2X0 = (Uf[r0] X0 - f0 (fp - Ff[r0]) dX0)/f0^2;
  scale = -c0/(4 w^2);
  chi = Xf[r] Dl/Sqrt[r^2 + a^2];
  Rup = 1/eta ((alpha + D[beta, r]/Dl) chi - beta/Dl D[chi, r]); dRup = D[Rup, r];
  rules = {Derivative[2][Xf][r] -> d2X0, Derivative[1][Xf][r] -> dX0, Xf[r] -> X0};
  scale {Rup /. rules /. r -> r0, dRup /. rules /. r -> r0}];
kerrPhotonMode[a_, l_, j_, nterms_: 120] := Module[{d = lrData[a, 30], r0, Ut, m, w, lam, S0, dS0, d2S0, Qf, Qin, k, rp, xA, rB = 40, rA, x0, YA, dYA, solIn, x, Y, r, rr, Yin0, dYin0, WY, pref, dpref, RIn, dRIn, RUp, dRUp, aIn, aUp, ZI, ZH},
  r0 = d["r0"]; Ut = d["Ups"]; m = l - j; w = N[m d["Omega"], 30];
  {lam, S0, dS0, d2S0} = swshData[a, l, m, w, nterms]; lam = N[lam]; w = N[w];
  Qf = QfunK[a, m, w, lam]; Qin = QwkbK[a, m, w, lam, -1];
  rp = rpK[a]; k = w - m a/(2 rp);
  (* horizon side: start at rA = r+ + 10^-6; Y_in -> (r+^2+a^2)^(1/2) Delta^(-s/2) e^(-i k rstar) times exp(Int g dr), with g finite at the horizon *)
  rA = rp + 10^-6; xA = N[rstarK[a, rA], 30];
  YA = N[Sqrt[rp^2 + a^2] Dl[a, rA]^(-s/2) Exp[-I k xA], 30] Exp[NIntegrate[(I (Qin[rr] + k) (rr^2 + a^2) + s (rr - 1))/Dl[a, rr], {rr, rp + 10^-12, rA}, PrecisionGoal -> 8, MaxRecursion -> 30]];
  dYA = I Qin[rA] YA; xA = N[xA];
  x0 = N[rstarK[a, r0]];
  solIn = NDSolve[{Y''[x] + Qf[r[x]] Y[x] == 0, r'[x] == drdrs[a, r[x]], Y[xA] == YA, Y'[xA] == dYA, r[xA] == N[rA]}, {Y, r}, {x, xA, x0},
     PrecisionGoal -> 12, AccuracyGoal -> 12, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9, "StiffnessTest" -> False}, InterpolationOrder -> All][[1]];
  {Yin0, dYin0} = {Y[x0], Y'[x0]} /. solIn;
  pref = Dl[a, r0]^(-s/2) (r0^2 + a^2)^(-1/2); dpref = D[Dl[a, rr]^(-s/2) (rr^2 + a^2)^(-1/2), rr] /. rr -> r0;
  RIn = pref Yin0; dRIn = dpref Yin0 + pref dYin0/drdrs[a, r0];
  {RUp, dRUp} = snUpAt[a, w, m, lam, N[r0], rB];
  WY = Dl[a, r0]^(s + 1) (RIn dRUp - RUp dRIn);   (* = Y_in Y_up' - Y_up Y_in' = 2 i w B_inc *)
  aIn = alphaSource[a, N[r0], 1, N[d["b"]], N[Ut], m, w, lam, N[S0], N[dS0], N[d2S0], RIn, dRIn];
  aUp = alphaSource[a, N[r0], 1, N[d["b"]], N[Ut], m, w, lam, N[S0], N[dS0], N[d2S0], RUp, dRUp];
  ZI = -8 Pi aIn/(WY Ut); ZH = -8 Pi aUp/(WY Ut);
  <|"a" -> a, "l" -> l, "j" -> j, "Flux" -> fluxKerr[a, m, w, lam, ZI, ZH]|>];
(* a = 0 must agree with photonMode; then the data of Fig. 2 at small l: mean flux summed over j <= 3 and +-m, compared with 2 g(a) Sum_{j<=3} kappa_j *)
Print["a = 0, l = 10, j = 0: Teukolsky ", kerrPhotonMode[0, 10, 0]["Flux"], "   Zerilli ", photonMode[10, 10]["Flux"]];
kerrList = {10, 20, 50};
kerrData = Table[Module[{res = Table[kerrPhotonMode[a, l, j], {j, 0, 3}]}, {a, l, l Total[(#["Flux"]["I"] + #["Flux"]["H"]) & /@ res]}], {a, {0, 1/2, -1/2, 9/10, -9/10}}, {l, kerrList}];
TableForm[Flatten[kerrData, 1], TableHeadings -> {None, {"a", "l", "l (Edot_I + Edot_H)"}}]
ListLogLinearPlot[Table[kerrData[[i, All, {2, 3}]], {i, 5}], Joined -> True, PlotMarkers -> Automatic, PlotLegends -> {"a=0", "a=0.5", "a=-0.5", "a=0.9", "a=-0.9"},
  Epilog -> Table[{Dashed, Line[{{Log[8], 2 Total[kappaKerr[a, #]["kappa"] & /@ Range[0, 3]]}, {Log[1200], 2 Total[kappaKerr[a, #]["kappa"] & /@ Range[0, 3]]}}]}, {a, {0, 1/2, -1/2, 9/10, -9/10}}],
  Frame -> True, FrameLabel -> {"l", "l Sum_m (Edot_I + Edot_H)/2"}] *)

(* ::Subsubsection:: *)
(*Kerr runs to l = 800 (commented out)*)

(* ::Text:: *)
(*This cell needs the previous one (and hence the Toolkit) and is commented out. Larger l, as in the paper: l = 100, 200, 400, 800 for the five spins, and the Richardson extrapolation in 1/l from l = 400 and 800 compared with the analytic 2 g(a) Sum_{j<=3} kappa_j (Sec. IV quotes agreement to better than 1.5 10^-5 in relative terms). Its output is the stored run read by the next cell.*)

(* kerrListLong = {100, 200, 400, 800};
kerrDataLong = Table[Module[{res = Table[kerrPhotonMode[a, l, j], {j, 0, 3}]}, {a, l, l Total[(#["Flux"]["I"] + #["Flux"]["H"]) & /@ res]}], {a, {0, 1/2, -1/2, 9/10, -9/10}}, {l, kerrListLong}];
TableForm[Flatten[kerrDataLong, 1], TableHeadings -> {None, {"a", "l", "l (Edot_I + Edot_H)"}}]
(* Richardson: with f(l) = f0 + c/l, f0 = 2 f(800) - f(400) *)
Table[With[{f400 = kerrDataLong[[i, 3, 3]], f800 = kerrDataLong[[i, 4, 3]], kA = 2 Total[kappaKerr[kerrDataLong[[i, 1, 1]], #]["kappa"] & /@ Range[0, 3]]},
   {kerrDataLong[[i, 1, 1]], 2 f800 - f400, kA, (2 f800 - f400)/kA - 1}], {i, 5}] // TableForm *)

(* ::Subsubsection:: *)
(*The stored Kerr runs: Richardson extrapolation and the last column of Table I*)

(* ::Text:: *)
(*The stored Kerr runs: kerr_mathematica_runL1.m and kerr_mathematica_runL2.m in notebook/ (the integrator of the previous cells, l = 100-800, rows {a, l, j, Edot_I, Edot_H}) and kerr_py/kerr_results.json (the Python implementation used for Fig. 2 and for the numbers of Sec. IV). Reproduces the Richardson extrapolation in 1/l from l = 400 and 800, which agrees with 2 g(a) Sum_{j<=3} kappa_j to better than 1.5 10^-5 in relative terms for the five spins (the omitted j >= 4 terms are 3 10^-5 of kappa, Section III C), the agreement between the two implementations, the last column of Table I (num./an. = 1.0012, 1.0014, 1.0017, 1.0011, 1.0010 at l = 800), and the O(l^(-1/2)) asymmetry between the two fluxes in Kerr, which was not computed analytically for a different from 0: Sqrt[l] (Edot_I - Edot_H)/(Edot_I + Edot_H) at l = 800 is 0.78, 1.10, 0.51 and 0.45 for a = 0.5, 0.9, -0.5 and -0.9, against 0.61 at a = 0 (the asymptotic sigma_bar of Section III C, which the a = 0 numerics reach by l = 800).*)

kerrMath = Join @@ (Get[FileNameJoin[{repoDir, "notebook", #}]] & /@ {"kerr_mathematica_runL1.m", "kerr_mathematica_runL2.m"});
kerrPy = Import[FileNameJoin[{repoDir, "kerr_py", "kerr_results.json"}], "JSON"];
mathSum[a_, l_] := With[{rows = Select[kerrMath, #[[1]] == a && #[[2]] == l &]}, l Total[rows[[All, 4]] + rows[[All, 5]]]];
pySum[a_, l_] := With[{rows = Select[kerrPy, ("a" /. #) == N[Abs[a]] && ("sign" /. #) == If[a >= 0, 1, -1] && ("l" /. #) == l && ("j" /. #) <= 3 &]}, l Total[(("FluxI" /. #) + ("FluxH" /. #)) & /@ rows]];
TableForm[Table[With[{f400 = mathSum[a, 400], f800 = mathSum[a, 800], kA = 2 Total[kappaKerr[a, #]["kappa"] & /@ Range[0, 3]]},
    {a, f400, f800, 2 f800 - f400, kA, (2 f800 - f400)/kA - 1, Max[Table[Abs[mathSum[a, l]/pySum[a, l] - 1], {l, {100, 200, 400, 800}}]]}], {a, {0, 1/2, -1/2, 9/10, -9/10}}],
  TableHeadings -> {None, {"a", "l=400", "l=800", "Richardson", "analytic", "rel. diff.", "max |Mathematica/Python - 1|"}}]
(* last column of Table I: numerical l (Edot_I + Edot_H) at l = 800, summed over m > 0 with j <= 3, over the analytic 2 g(a) Sum_{j<=3} kappa_j *)
Print["Table I, num./an. at l = 800: ", Table[{a, N[pySum[a, 800]/(2 Total[kappaKerr[a, #]["kappa"] & /@ Range[0, 3]]), 5]}, {a, {0, 1/2, 9/10, -1/2, -9/10}}]];
(* the O(l^(-1/2)) asymmetry between the two fluxes in Kerr, from the Python data: Sqrt[l] (I - H)/(I + H) tends to a spin-dependent constant *)
pyIH[a_, l_] := With[{rows = Select[kerrPy, ("a" /. #) == N[Abs[a]] && ("sign" /. #) == If[a >= 0, 1, -1] && ("l" /. #) == l && ("j" /. #) <= 3 &]}, {Total[("FluxI" /. #) & /@ rows], Total[("FluxH" /. #) & /@ rows]}];
TableForm[Table[Prepend[Table[With[{ih = pyIH[a, l]}, Sqrt[l] (ih[[1]] - ih[[2]])/(ih[[1]] + ih[[2]])], {l, {100, 200, 400, 800}}], a], {a, {0, 1/2, 9/10, -1/2, -9/10}}], TableHeadings -> {None, {"a", "l=100", "l=200", "l=400", "l=800"}}]
Print["Sqrt[l] (I - H)/(I + H) at l = 800 for a = 0.5, 0.9, -0.5, -0.9: ", Table[With[{ih = pyIH[a, 800]}, N[Sqrt[800] (ih[[1]] - ih[[2]])/(ih[[1]] + ih[[2]]), 3]], {a, {1/2, 9/10, -1/2, -9/10}}], "  (paper: 0.78, 1.10, 0.51, 0.45);  a = 0: ", With[{ih = pyIH[0, 800]}, N[Sqrt[800] (ih[[1]] - ih[[2]])/(ih[[1]] + ih[[2]]), 3]], " against the asymptotic sigma_bar = ", N[sigBar, 3]];

(* ::Subsubsection:: *)
(*Fig. 2*)

(* ::Text:: *)
(*Fig. 2 from the stored data: top, l (Edot_I + Edot_H)/2 summed over the four dominant m (j <= 3) and over +-m from kerr_results.json, with the analytic 2 g(a) Sum_{j<=3} kappa_j dashed; bottom, the spin factor g(a) with the five spins marked.*)

fig2spins = {9/10, 1/2, 0, -1/2, -9/10};
fig2top = Table[Table[{l, pySum[a, l]}, {l, {10, 20, 50, 100, 200, 400, 800}}], {a, fig2spins}];
fig2lines = Table[2 Total[kappaKerr[a, #]["kappa"] & /@ Range[0, 3]], {a, fig2spins}];
GraphicsColumn[{
  Show[ListLogLinearPlot[fig2top, Joined -> True, PlotMarkers -> Automatic, PlotLegends -> ("a/M = " <> ToString[N[#]] & /@ fig2spins), PlotRange -> {{8, 1200}, {0.04, 0.1}}, Frame -> True, FrameLabel -> {"l", "l (Edot_I + Edot_H) M^2/E^2"}],
    Graphics[{Dashed, Line[{{Log[8], #}, {Log[1200], #}}] & /@ fig2lines}]],
  ListLinePlot[gCurve, Frame -> True, FrameLabel -> {"a/M", "g(a)"}, PlotRange -> {{-1, 1}, {0.5, 1.5}}, GridLines -> {None, {1}},
    Epilog -> {PointSize[Medium], Point[Table[{N[a], kappaKerr[a, 0]["kappa"]/kappa00}, {a, fig2spins}]]}]}]

(* ::Section:: *)
(*V. Comparison with the geodesic synchrotron radiation literature (and the Chrzanowski-Misner comparison of Sec. IV)*)

(* ::Subsubsection:: *)
(*Eqs. (22)-(23) and their null limits*)

(* ::Text:: *)
(*Eqs. (22)-(23): the formulae of Breuer, Ruffini, Tiomno and Vishveshwara (PRD 7, 1002, 1973, Eqs. 14-15) for a massive particle at r0 = 3 (1 + delta), per unit mu^2 and in the notation of the paper (eta = 1/2 + 3 m delta/2, eta' = eta + 1). Their null limit, mu^2/delta -> 9 E^2 at fixed E, gives exactly kappa_0 for the even mode and exactly kappa_1/4 for the odd one; the cell takes the limit in closed form and checks it numerically.*)

Peven[m_, delta_] := With[{eta = 1/2 + 3 m delta/2}, Exp[-Pi eta/2]/(162 Pi^(3/2) m delta) Abs[(eta + 1/2) Gamma[1/4 + I eta/2] + Sqrt[2] (1 - I)/Sqrt[3 m] Gamma[3/4 + I eta/2]]^2];   (* Eq. (22) *)
Podd[m_, delta_] := With[{eta = 3/2 + 3 m delta/2}, Exp[-Pi eta/2]/(162 Pi^(3/2) m delta) Abs[Sqrt[2] Gamma[3/4 + I eta/2] + (1 + I)/(2 Sqrt[3 m]) Gamma[1/4 + I eta/2]]^2];           (* Eq. (23) *)
(* null limit: multiply by mu^2/E^2 = 9 delta, delta -> 0, m -> Infinity, times m *)
(* closed-form limits (the cross terms vanish as m^(-1/2)), and a direct numerical evaluation at m delta = 10^-8, m = 10^12 *)
limEven = 9 Exp[-Pi/4]/(162 Pi^(3/2)) Abs[Gamma[1/4 + I/4]]^2;
limOdd = 9 Exp[-3 Pi/4]/(162 Pi^(3/2)) 2 Abs[Gamma[3/4 + 3 I/4]]^2;
Print["even: limit/kappa_0 = ", N[limEven/kappaEven[0], 15], ",  numerical: ", N[9 delta m Peven[m, delta]/kappaEven[0] /. {delta -> 10^-20, m -> 10^12}, 15]];
Print["odd:  limit/kappa_1 = ", N[limOdd/kappaOdd[1], 15], ",  numerical: ", N[9 delta m Podd[m, delta]/kappaOdd[1] /. {delta -> 10^-20, m -> 10^12}, 15], "   (= 1/4)"];

(* ::Subsubsection:: *)
(*Test of the odd-parity formula at finite m*)

(* ::Text:: *)
(*Test of Eqs. (22)-(23) at finite m, away from the m -> Infinity limit (Sec. V): exact fluxes to infinity for a timelike circular orbit at r0 = 3 (1 + delta), delta = 10^-4, for the even mode l = m and the odd mode l = m + 1, m = 40, 100, 200, 400 (per unit mu^2). The paper quotes Teukolsky values; here the Zerilli and Regge-Wheeler equations with the massive source are integrated with the same algorithm as in Section III D (circMode), which gives the same fluxes. The even formula is accurate to O(m^(-1/2)); the ratio of the exact odd flux to Eq. (23) is 3.15, 3.48, 3.64, 3.77, and fits rho_inf + rho_1 m^(-1/2) + rho_2/m and rho_inf + rho_1 m^(-1/2) to the four values give rho_inf = 4.06 and 4.05: the ratio tends to 4 with O(m^(-1/2)) corrections. With the 1973 odd term, kappa would be 0.0574 (kappa_1 replaced by kappa_1/4, j <= 1 only) or 0.0580 (all odd kappa_j divided by four), against 0.0637 and the numerical fit 0.064 +- 0.001 of the 2021 paper.*)

circMode[l_, m_, r0_, rsA_: -40, rB_: 40, pg_: 12] := Module[{parity, w, L0, E0, V, Qin, Qup, rA, ampA, ampB, psiA, dpsiA, psiB, dpsiB, rsB, rs0, eqs, solIn, solUp, r, Psi, rr, PsiIn, dPsiIn, PsiUp, dPsiUp, W, J, P, dP, ZI, ZH},
  parity = If[EvenQ[l + m], "Zerilli", "ReggeWheeler"];
  w = N[m Sqrt[1/r0^3]]; L0 = N[Sqrt[r0^2/(r0 - 3)], 30]; E0 = N[Sqrt[(r0 - 2) (r0^2 + L0^2)/r0^3], 30];   (* E0 = E/mu = gamma, L0 = L/mu *)
  V = Vpot[l, parity, #] &;
  rA = N[r /. FindRoot[rstar[r] == rsA, {r, 2 + 2 Exp[(rsA - 2)/2], 2 + 10^-30, 3}, WorkingPrecision -> 40, AccuracyGoal -> 35], 40];
  Qin = QWKB[l, parity, w, -1]; Qup = QWKB[l, parity, w, +1];
  ampA = Exp[-NIntegrate[Im[Qin[rr]]/f[rr], {rr, 2, rA}, PrecisionGoal -> 12, MaxRecursion -> 30]];
  ampB = Exp[NIntegrate[Im[Qup[rr]]/f[rr], {rr, rB, Infinity}, PrecisionGoal -> 12, MaxRecursion -> 30]];
  psiA = ampA; dpsiA = I Qin[rA] ampA; psiB = ampB; dpsiB = I Qup[rB] ampB;
  rsB = rstar[rB]; rs0 = rstar[r0];
  eqs = {Psi''[rr] + (w^2 - V[r[rr]]) Psi[rr] == 0, r'[rr] == f[r[rr]]};
  solIn = NDSolve[Join[eqs, {Psi[rsA] == psiA, Psi'[rsA] == dpsiA, r[rsA] == rA}], {Psi, r}, {rr, rsA, rs0},
     PrecisionGoal -> pg, AccuracyGoal -> pg, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}, InterpolationOrder -> All][[1]];
  solUp = NDSolve[Join[eqs, {Psi[rsB] == psiB, Psi'[rsB] == dpsiB, r[rsB] == rB}], {Psi, r}, {rr, rs0, rsB},
     PrecisionGoal -> pg, AccuracyGoal -> pg, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}, InterpolationOrder -> All][[1]];
  PsiIn = Psi[rs0] /. solIn; dPsiIn = Psi'[rs0] /. solIn; PsiUp = Psi[rs0] /. solUp; dPsiUp = Psi'[rs0] /. solUp;
  W = PsiIn dPsiUp - PsiUp dPsiIn;
  J = jumpsZRW[l, m, r0, E0, L0, Conjugate[YLRnum[l, m]], Conjugate[dYLRnum[l, m]]]; P = N[J[[1]]]; dP = N[f[r0] J[[2]]];
  ZI = (PsiIn dP - P dPsiIn)/W; ZH = (PsiUp dP - P dPsiUp)/W;
  <|"l" -> l, "m" -> m, "r0" -> r0, "gamma" -> E0, "Flux" -> fluxFromZ[l, m, w, ZI, ZH]|>];
delta4 = 10^-4;
brtvTest = Table[{m, circMode[m, m, 3 (1 + delta4)]["Flux"]["I"]/Peven[m, delta4], circMode[m + 1, m, 3 (1 + delta4)]["Flux"]["I"]/Podd[m, delta4]}, {m, {40, 100, 200, 400}}];
TableForm[brtvTest, TableHeadings -> {None, {"m", "even: exact/1973", "odd: exact/1973"}}]
(* extrapolation of the odd ratio in powers of m^(-1/2) *)
fitOdd3 = Fit[brtvTest[[All, {1, 3}]], {1, 1/Sqrt[m], 1/m}, m]; fitOdd2 = Fit[brtvTest[[All, {1, 3}]], {1, 1/Sqrt[m]}, m];
Print["odd ratio, m -> Infinity: three-term fit ", fitOdd3 /. m -> Infinity, ",  two-term fit ", fitOdd2 /. m -> Infinity, "  (paper: 4.06, 4.05)"];
Print["kappa with the 1973 odd term: j <= 1 only: ", N[2 (kappaEven[0] + kappaOdd[1]/4), 5], ";  all odd j divided by 4: ", N[2 Sum[If[EvenQ[j], kappaJ[j], kappaJ[j]/4], {j, 0, 20}], 5], "  (paper: 0.0574, 0.0580)"];

(* ::Subsubsection:: *)
(*Chrzanowski-Misner: the null limit*)

(* ::Text:: *)
(*Kerr, the paragraph after Eq. (21): the gravitational power formula of Chrzanowski and Misner (PRD 10, 1701, 1974), their Eq. (4.28), with their Eq. (2.35b) for the cut-off harmonic m_crit and Eq. (2.31) for the energy per unit rest mass gamma = E/mu, in the null limit gamma -> infinity at fixed E (their epsilon -> 1, k = 0). The power per mode is (12/Sqrt[Pi]) E^(-Pi/2) E^2 (r0 - M)/[r0^2 (r0 + 3M) m]: its spin dependence coincides exactly with the closed form of g(a), Eq. (21), and its constant is 1.877 times kappa_0, i.e. 0.94 times 2 kappa_0, the flux per |m| (both signs of m), which is the quantity their formula refers to (next cell).*)

mcritCM[r_, gam2_] := (2 Sqrt[3]/Pi) (r + 3)/Sqrt[r] gam2;                                           (* CM Eq. (2.35b), M = 1 *)
PCM[r_, m_, gam2_, eps_] := 2 Sqrt[Pi] (1/gam2) (r - 1) Sqrt[3 r]/(r^2 (r + 3)^2) eps^(3/2) (mcritCM[r, gam2]/m) Exp[-Pi eps/2];   (* CM Eq. (4.28), k = 0 term, per E^2 = mu^2 gamma^2 *)
nullCM[r_] := Simplify[m PCM[r, m, g2, 1]];
Print["Chrzanowski-Misner, null limit: m P_m/E^2 = ", nullCM[r], ";  ratio to (r0-1)/(r0^2 (r0+3)): ", Simplify[nullCM[r]/((r - 1)/(r^2 (r + 3)))], " = (12/Sqrt[Pi]) E^(-Pi/2) = ", N[12/Sqrt[Pi] Exp[-Pi/2], 6]];
Print["spin factor nullCM[r0]/nullCM[3] - g(a) of Eq. (21): ", Simplify[nullCM[r]/nullCM[3] - 27 (r - 1)/(r^2 (r + 3))]];
Print["Schwarzschild constant: CM ", N[nullCM[3], 8], " vs kappa_0 = ", N[kappaJ[0], 8], ",  ratio = ", N[nullCM[3]/kappaJ[0], 6]];

(* ::Subsubsection:: *)
(*The m-counting of the 1973 formulae: a scalar-flux check*)

(* ::Text:: *)
(*Which quantity the 1973-74 formulae refer to. The power 'at frequency m omega_0' of Breuer, Chrzanowski, Hughes and Misner (PRD 8, 4309, 1973, Eq. 5.4) and of Chrzanowski and Misner collects the modes m and -m. This cell checks it directly: the exact flux at infinity of a scalar charge q (field equation Box Phi = 4 pi q Int dtau delta^4/Sqrt[-g], stress tensor (1/4 pi)(dPhi dPhi - g (dPhi)^2/2), as in BCHM) on a circular orbit at r0 = 3 + delta', mode l = m, single m > 0, computed with the scalar Regge-Wheeler equation u'' + (omega^2 - V0) u = S0 delta(r* - r0* ), V0 = f (l(l+1)/r^2 + 2/r^3), S0 = 4 pi Ybar_lm(pi/2, 0)/(u^t r0), flux (1/4 pi) omega^2 |Z|^2. BCHM Eq. (5.4) is twice the single-m flux (the residual deviation is the O(m^(-1/2)) error of their asymptotic formula).*)

scalarFlux[l_, r0_, rB_: 2000] := Module[{Om, w, ut, m = l, S0, V, eqs, rA, rsA, solIn, rsB, solUp, uin, duin, uup, duup, W, Z, rs0},
  Om = Sqrt[1/r0^3]; w = m Om; ut = 1/Sqrt[1 - 3/r0]; S0 = 4 Pi Conjugate[SphericalHarmonicY[l, m, Pi/2, 0]]/(ut r0);
  V[r_] := f[r] (l (l + 1)/r^2 + 2/r^3);
  eqs = {u''[x] == -(w^2 - V[r[x]]) u[x], r'[x] == f[r[x]]};
  rsA = -40; rA = r /. FindRoot[rstar[r] == rsA, {r, 2 + 2 Exp[(rsA - 2)/2], 2 + 10^-30, 3}, WorkingPrecision -> 40];
  solIn = NDSolve[Join[eqs, {u[rsA] == Exp[-I w rsA], u'[rsA] == -I w Exp[-I w rsA], r[rsA] == rA}], {u, r}, {x, rsA, rstar[r0]}, PrecisionGoal -> 12, AccuracyGoal -> 12, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}][[1]];
  rsB = rstar[rB];   (* outgoing data imposed where V/omega^2 is negligible *)
  solUp = NDSolve[Join[eqs, {u[rsB] == Exp[I w rsB], u'[rsB] == I w Exp[I w rsB], r[rsB] == rB}], {u, r}, {x, rstar[r0], rsB}, PrecisionGoal -> 12, AccuracyGoal -> 12, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}][[1]];
  rs0 = rstar[r0]; {uin, duin} = {u[rs0], u'[rs0]} /. solIn; {uup, duup} = {u[rs0], u'[rs0]} /. solUp;
  W = uin duup - duin uup; Z = uin S0/W; w^2 Abs[Z]^2/(4 Pi)];
bchm54[m_, deltaB_] := With[{mcrit = 4/(Pi deltaB)}, With[{eps = 1 + (4/Pi) m/mcrit}, (1/(27 Pi^(5/2))) (m/mcrit) Exp[-Pi eps/4] Abs[Gamma[1/4 + I eps/4]]^2]];
(* normalization check in the weak field: a scalar charge on a circular orbit of radius r0 >> M radiates the dipole power q^2 a^2/3 = q^2 M^2/(3 r0^4), shared equally by m = +1 and m = -1, so the single-m flux must tend to M^2/(6 r0^4) *)
Print["weak-field check, l = m = 1: single-m flux / [M^2/(6 r0^4)] = ", Table[{r0, scalarFlux[1, r0, 2 10^6] 6 r0^4}, {r0, {100, 300}}], "  (-> 1 with O(M/r0) corrections)"];
TableForm[Table[With[{num = scalarFlux[l, 3 + dB], an = bchm54[l, dB]}, {l, N[dB], N[1 + l dB], num, N[an], N[an/num]}], {l, {40, 100}}, {dB, {1/100}}] // Flatten[#, 1] &,
  TableHeadings -> {None, {"l = m", "delta'", "epsilon", "single-m scalar flux", "BCHM (5.4)", "BCHM/(single m)"}}]

(* ::Subsubsection:: *)
(*Origin of the 0.94*)

(* ::Text:: *)
(*Origin of the 0.94 (same paragraph of the paper): two errors of the WKB treatment compensate. In the massive regime (large eta) the Chrzanowski-Misner master formula, their Eq. (5.1) with (5.2), agrees with the exact scalar result of BCHM, Eq. (5.3), but tends to 8 times the exact tensor result of Breuer, Ruffini, Tiomno and Vishveshwara, Eq. (22), which is per single m: referred to the same counting of m, their tensor formula is 4 times the exact one, the (s!)^2 of their master formula. In the null limit this is compensated by the WKB form of the barrier-top factor, eta^(3/2) E^(-Pi eta), which is the large-eta form of the exact (eta + 1/2)^2 F(eta) obtained with Stirling's approximation of the Gamma function (F(eta) -> E^(-Pi eta)/(Pi Sqrt[eta/2]), Sec. VI A): at eta = 1/2 it underestimates the exact factor by 4 x 1.065, i.e. the ratio WKB/exact is 0.2347 = 1/(4 x 1.065). What remains per |m| is 4 x 0.2347 = 1/1.065 = 0.94 (and 8 x 0.2347 = 1.877 of the single-m kappa_0). This factor 4 is unrelated to the odd-parity factor 4 of Sec. V. Conventions: Schwarzschild, M = 1, per unit mu^2, l = m; the delta of BCHM and BRTV is 3 times ours, gamma^2 = 1/(9 delta), m_crit = 4/(Pi delta') = 12 gamma^2/Pi, epsilon = 1 + 3 m delta = 2 eta.*)

epsCM[m_, d_] := 1 + 3 m d; mcritCM0[d_] := 4/(3 Pi d);
PBCHM[m_, d_] := (1/(27 Pi^(5/2))) (m/mcritCM0[d]) Exp[-Pi epsCM[m, d]/4] Abs[Gamma[1/4 + I epsCM[m, d]/4]]^2;     (* BCHM Eq. (5.3), exact scalar, q = 0 *)
ACM3 = 2 (3 - 1) Sqrt[3 3]/(Sqrt[Pi] 3^2 (3 + 3)^2);                                                             (* CM Eq. (5.2) at r_gamma = 3 *)
PCMs[s_, m_, d_] := ACM3 (s!)^2 (4 m/(epsCM[m, d] Pi mcritCM0[d]))^(1 - s) Sqrt[epsCM[m, d]] Exp[-Pi epsCM[m, d]/2];   (* CM Eq. (5.1), k = 0 term *)
Print["CM Eq. (4.28), k = 0, equals the master formula with s = 2: ", Simplify[PCMs[2, m, d]/(2 Sqrt[Pi] (3 - 1) Sqrt[3 3]/(3^2 (3 + 3)^2) epsCM[m, d]^(3/2) (mcritCM0[d]/m) Exp[-Pi epsCM[m, d]/2])]];
TableForm[Table[With[{d = 10^-3}, {m, N[epsCM[m, d]/2], N[PCMs[0, m, d]/PBCHM[m, d], 6], N[PCMs[2, m, d]/Peven[m, d], 6], N[Peven[m, d]/PBCHM[m, d], 6], N[PCMs[2, m, d]/PCMs[0, m, d], 6]}], {m, {100, 1000, 10000, 100000, 1000000}}],
  TableHeadings -> {None, {"m", "eta", "CM scalar/BCHM", "CM tensor/BRTV", "BRTV/BCHM (exact)", "CM tensor/scalar"}}]
Print["null limit: CM scalar/BCHM = ", N[Limit[PCMs[0, m, d]/PBCHM[m, d] /. d -> 0, m -> Infinity], 6], ",  CM tensor/BRTV = ", N[Limit[PCMs[2, m, d]/Peven[m, d] /. d -> 0, m -> Infinity], 6]];
wkbOverExact = N[(Sqrt[2] (1/2)^(3/2) Exp[-Pi/2]/Pi)/((1/2 + 1/2)^2 Fpar[1/2]), 6];   (* [Sqrt[2] eta^(3/2) E^(-Pi eta)/Pi] / [(eta+1/2)^2 F(eta)] at eta = 1/2 *)
Print["WKB/exact barrier factor at eta = 1/2: ", wkbOverExact, " = 1/", 1/wkbOverExact, " = 1/(4 x ", 1/(4 wkbOverExact), ");  per |m|: 4 x ", wkbOverExact, " = ", 4 wkbOverExact, " = 1/1.065 (paper: 0.94);  per single m: 8 x ", wkbOverExact, " = ", 8 wkbOverExact];

(* ::Section:: *)
(*VI. Discussion*)

(* ::Subsection:: *)
(*A. Timelike orbits near the light ring*)

(* ::Subsubsection:: *)
(*Ingredients of Eq. (24)*)

(* ::Text:: *)
(*The ingredients of Eq. (24): for a particle at r0 = 3 (1 + delta), gamma^2 = E^2/mu^2 = (1 + 3 delta)^2/[9 delta (1 + delta)], and the large-l expansion with delta = sigma/l gives, in place of Eq. (5), [dPsi/dr]/E = -8 Pi Y (1 + 3 l delta/2)/(l M) and the barrier parameter eta = 1/2 + 3 l delta/2 for m = l.*)

Clear[l, sig, dd, r];
r0t = 3 (1 + dd); L0t = Sqrt[r0t^2/(r0t - 3)]; E0t = Sqrt[(r0t - 2) (r0t^2 + L0t^2)/r0t^3];
Print["gamma^2 = ", Simplify[E0t^2]];
dPsit = jumpsZRW[l, l, r0t, E0t, L0t, Y, 0, True][[2]];
Print["[Psi']/E (m = l, delta = sigma/l) = ", Series[(dPsit/E0t) /. dd -> sig/l, {l, Infinity, 1}] // Normal // Simplify];
Print["eta (m = l, delta = sigma/l) = ", Series[(VZ[l, 3] - l^2/r0t^3)/Sqrt[2 kZ] /. dd -> sig/l, {l, Infinity, 0}] // Normal // Simplify, "   (omega^2 = l^2/r0^3 for the massive orbit, V0 and k as in Section III B)"];

(* ::Subsubsection:: *)
(*Eq. (24) against the timelike fluxes*)

(* ::Text:: *)
(*Eq. (24), Edot_ll = (kappa_0/l) (eta + 1/2)^2 F(eta)/F(1/2), against the numerical fluxes of timelike orbits (Sec. VI A). The stored file timelike_results.m in the repository root contains the fluxes of the dominant mode l = m for orbits with delta = 1/(9 gamma^2), gamma = 5, 10, 20, computed with the integrator of Section III D for the massive source (rows {r0, l, Edot_I/E^2, Edot_H/E^2}); since Eq. (24) is the leading term in l^(-1/2), it is compared with the mean of the two fluxes, in which the O(l^(-1/2)) correction cancels. As stated in Sec. VI A, the agreement is better than 1% for gamma = 20 and 120 <= l <= 4800, 0.4-4% for gamma = 10 and 30 <= l <= 1200, and better than 4% for gamma = 5 and 25 <= l <= 125 (the table lists all the stored rows, including the smaller l, where the O(1/l) corrections are larger, and the larger l of the gamma = 5 and 10 runs, deep in the exponential cut-off); the gamma column gives the nominal gamma = 1/Sqrt[9 delta] of the runs. A few of the rows are then recomputed directly with circMode. Then the cut-off: F(eta) ~ E^(-Pi eta)/(Pi Sqrt[eta/2]) at large eta, so that with 3 l delta/2 = l/(6 gamma^2) the flux is cut off as (l/(6 gamma^2))^(3/2) Exp[-Pi l/(6 gamma^2)], i.e. c_1 = Pi/6 = 0.52 in the notation of the 2021 paper (a pure exponential fitted over a finite range of l returns a smaller effective coefficient, 0.42 +- 0.02 there); and the sum over l up to the cut-off l ~ gamma^2, which gives for each flux Edot_tot M^2/E^2 = kappa Log[gamma^2] + const = 2 kappa Log[gamma] = 0.127 Log[gamma], consistent with the slope k_1 = 0.12 +- 0.01 fitted in the 2021 paper.*)

Edotll[l_, delta_] := With[{eta = 1/2 + 3 l delta/2}, kappaEven[0]/l (eta + 1/2)^2 Fpar[eta]/Fpar[1/2]];   (* Eq. (24) *)
gammaOf[delta_] := Sqrt[(1 + 3 delta)^2/(9 delta (1 + delta))];
timelikeRun = Get[FileNameJoin[{repoDir, "timelike_results.m"}]];
TableForm[Map[With[{d = #[[1]]/3 - 1, l = #[[2]], mean = (#[[3]] + #[[4]])/2}, {Round[1/Sqrt[9 d]], Round[gammaOf[d], 0.01], l, 1/2 + 3 l d/2, l mean, l Edotll[l, d], mean/Edotll[l, d] - 1}] &, timelikeRun],
  TableHeadings -> {None, {"gamma (nominal)", "gamma (exact)", "l", "eta", "l (Edot_I + Edot_H)/2 / E^2", "l Eq. (24)", "rel. diff."}}]
(* the ranges quoted in Sec. VI A: smallest and largest relative deviation of the mean flux from Eq. (24) *)
rangeCheck[gamma_, lmin_, lmax_] := MinMax[Abs[(#[[3]] + #[[4]])/2/Edotll[#[[2]], #[[1]]/3 - 1] - 1] & /@ Select[timelikeRun, Round[1/Sqrt[9 (#[[1]]/3 - 1)]] == gamma && lmin <= #[[2]] <= lmax &]];
Print["relative deviation from Eq. (24), {min, max}: gamma = 20, 120 <= l <= 4800: ", rangeCheck[20, 120, 4800], " (paper: better than 1%);  gamma = 10, 30 <= l <= 1200: ", rangeCheck[10, 30, 1200], " (paper: 0.4-4%);  gamma = 5, 25 <= l <= 125: ", rangeCheck[5, 25, 125], " (paper: better than 4%)"];
(* direct recomputation of a few rows (per unit mu^2 from circMode, divided by gamma^2 to express them per unit E^2) *)
TableForm[Table[With[{gamma = gl[[1]], l = gl[[2]]}, With[{d = 1/(9 gamma^2), md = circMode[gl[[2]], gl[[2]], 3 (1 + 1/(9 gamma^2))]},
    With[{num = (md["Flux"]["I"] + md["Flux"]["H"])/2/md["gamma"]^2, stored = First[Select[timelikeRun, Abs[#[[1]] - 3 (1 + 1/(9 gamma^2))] < 10^-6 && #[[2]] == l &]]},
      {gamma, l, num, num/((stored[[3]] + stored[[4]])/2) - 1, num/Edotll[l, d] - 1}]]], {gl, {{5, 50}, {5, 125}, {10, 100}, {10, 300}, {20, 120}, {20, 400}}}],
  TableHeadings -> {None, {"gamma", "l", "numerical (mean), per E^2", "rel. diff. from the stored run", "rel. diff. from Eq. (24)"}}]
Print["F(eta)/(E^(-Pi eta)/(Pi Sqrt[eta/2])) at eta = 10, 100, 1000: ", Table[N[Fpar[eta]/(Exp[-Pi eta]/(Pi Sqrt[eta/2])) /. eta -> 10^k, 10], {k, 1, 3}]];
Print["l Edot_ll/kappa_0 divided by x^(3/2) Exp[-Pi x], x = l/(6 gamma^2) = eta - 1/2, at x = 5, 10, 20, 40 (tends to a constant): ", Table[N[(x + 1)^2 Fpar[x + 1/2]/Fpar[1/2]/(x^(3/2) Exp[-Pi x]), 8], {x, {5, 10, 20, 40}}]];
Print["cut-off coefficient c_1 = Pi/6 = ", N[Pi/6, 4], "  (pure-exponential fit of the 2021 paper: 0.42 +- 0.02)"];
Print["coefficient of Log[gamma] in Edot_tot M^2/E^2 (each flux): kappa Log[gamma^2] = 2 kappa Log[gamma] = 4 Sum_j kappa_j Log[gamma] = ", N[2 kappaTot, 4], " Log[gamma]  (fit of the 2021 paper: k_1 = 0.12 +- 0.01)"];

(* ::Section:: *)
(*VII. Conclusions: energy per orbit and the M87* estimate*)

(* ::Subsubsection:: *)
(*Energy per orbit and the M87* estimate*)

(* ::Text:: *)
(*Numbers of Sec. VII: multiplied by the orbital period 2 Pi b = 6 Sqrt[3] Pi M and summed over the two fluxes, Eq. (13) gives the energy radiated per orbit by an ultrarelativistic body near the critical impact parameter, Delta E = 2 kappa (2 Pi b) (E^2/M) Log[l_max] = 4.2 (E^2/M) Log[l_max], half of it absorbed, in the regime Delta E << E. Then the estimate for the millimetre photons of the Event Horizon Telescope around M87* (M = 6.5 10^9 solar masses, wavelength lambda = 1.3 mm): the wavelength cuts the sum at l_max ~ E r0/hbar = 2 Pi r0/lambda ~ 10^17, E/M ~ 10^-79, and the fractional energy loss per orbit is Delta E/E = 4.2 (E/M) Log[l_max] ~ 10^-77. The last estimate, for Sgr A*, is not quoted in the paper.*)

Print["orbital period 2 Pi b = 6 Sqrt[3] Pi = ", N[2 Pi bLR, 4], ";  Delta E per orbit = 2 kappa (2 Pi b) E^2/M = ", N[2 kappaTot 2 Pi bLR, 4], " (E^2/M) Log[l_max]   (both fluxes)"];
With[{G = 6.674 10^-8, c = 2.998 10^10, h = 6.626 10^-27, Msun = 1.989 10^33, lambda = 0.13},
  Module[{M = 6.5 10^9 Msun, EoverM, lmax},
    EoverM = (h c/lambda)/(M c^2);               (* photon energy over the rest energy of the hole = E/M in geometric units *)
    lmax = 2 Pi (3 G M/c^2)/lambda;              (* E r0/hbar = 2 Pi r0/lambda, r0 = 3 M *)
    Print["M87*: E/M = ", EoverM, ",  l_max = 2 Pi r0/lambda = ", lmax, ",  Delta E/E per orbit = ", N[2 kappaTot 2 Pi bLR] EoverM Log[lmax]]]];
(* the same estimate in physical units for Sgr A* (4 10^6 solar masses, 1 mm photons, a luminosity of 10^36 erg/s), not quoted in the paper: GW power per photon and incoherent total *)
With[{G = 6.674 10^-8, c = 2.998 10^10, h = 6.626 10^-27, Msun = 1.989 10^33},
  Module[{M = 4.0 10^6 Msun G/c^2, Ephot = h c/0.1, EG, lmax, Nph, Pph},
    EG = G Ephot/c^4;                         (* photon energy in cm *)
    lmax = 2 Pi 3 M/0.1;                      (* E r0/hbar = 2 Pi r0/wavelength *)
    Pph = 2 kappaTot EG^2/M^2 Log[lmax] c^5/G;  (* erg/s per photon, both fluxes *)
    Nph = 10^36/Ephot (2 Pi 3 Sqrt[3] M/c);   (* photons on the ring during one orbit at 10^36 erg/s *)
    Print["Sgr A*: per photon ", Pph, " erg/s;  N_photons ~ ", Nph, ";  incoherent total ~ ", Nph Pph, " erg/s"]]];

(* ::Section:: *)
(*Summary of the numbers quoted in the paper*)

(* ::Text:: *)
(*Sec. III B: V0 = l(l+1)/27, k = 2 l^2/729, V3 = 4 l^2/6561; eta_j = j + 1/2; lambda_L = Sqrt[2k]/(2 omega) = Omega; |T|^2 = 0.041 at eta = 1/2 (96% reflected). Sec. III C: kappa_0 = 0.0277668, kappa_1 = 0.0037670, kappa_2 = 2.7521 10^-4, kappa_3 = 1.6604 10^-5, kappa_4 = 9.23 10^-7; Sum_j kappa_j = 0.0318265; kappa = 0.063653; cross term 1.039, cubic term -0.23, sigma_0 = 0.81, sigma_j < 0 for j >= 1, sigma_bar = 0.61; 2 sigma_bar/Sqrt[l] = 12% at l = 100; horizon fraction 40% for l_max = 100-140. Sec. III D: j <= 4 captures all but 2 10^-6 of the sum; three-term fit 0.0318263 against 0.0318265; at l = 12800 the kappa_j, j <= 4, to better than 10^-3, kappa_0 to 10^-4 and kappa_1 to 2 10^-5; at l = 6400 the cross term 1.0389 (analytic 1.0390) and the barrier asymmetry -0.231. Sec. IV and Table I: r0, b, Upsilon_t, k~, |Ahat_0|, N, kappa_0(a) and g(a) for a = 0, 0.5, 0.9, 0.99, -0.5, -0.9, -0.99; N = 16/r0^6; g(a) = 27 (r0 - 1)/(r0^2 (r0 + 3)), maximum 6 Sqrt[3] - 9 = 1.392 at r0 = Sqrt[3], a = 0.834, g = 81/112 = 0.723 at r0 = 4, g = 0.99 at a = 0.98 and 0.62 at a = 0.995; horizon/infinity ratio 1; num./an. = 1.0012, 1.0014, 1.0017, 1.0011, 1.0010 at l = 800; Richardson extrapolation within 1.5 10^-5 (omitted j >= 4 terms: 3 10^-5 of kappa); Kerr asymmetries 0.78, 1.10, 0.51, 0.45 at l = 800 against 0.61 at a = 0. Chrzanowski-Misner: (12/Sqrt[Pi]) E^(-Pi/2) E^2 (r0 - 1)/(r0^2 (r0 + 3) m), constant 0.94 of 2 kappa_0 (1.877 kappa_0); their tensor formula (s!)^2 = 4 times the exact one in the massive regime, compensated at eta = 1/2 by the WKB barrier factor, 1/(4 x 1.065) = 0.2347. Sec. V: null limit of Eqs. (22)-(23) = kappa_0 and kappa_1/4; odd-parity ratio exact/1973 = 3.15, 3.48, 3.64, 3.77 for m = 40, 100, 200, 400 at delta = 10^-4, fits 4.06 and 4.05; kappa with the 1973 odd term 0.0574 (j <= 1) and 0.0580 (all odd j), against 0.0637. Sec. VI A: Eq. (24) within 1% (gamma = 20, 120 <= l <= 4800), 0.4-4% (gamma = 10, 30 <= l <= 1200), better than 4% (gamma = 5, 25 <= l <= 125); c_1 = Pi/6; Edot_tot M^2/E^2 = 0.127 Log[gamma]. Sec. VII: Delta E per orbit = 4.2 (E^2/M) Log[l_max]; M87*: l_max ~ 10^17, E/M ~ 10^-79, Delta E/E ~ 10^-77 per orbit.*)
