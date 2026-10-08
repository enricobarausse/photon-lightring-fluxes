(* ::Package:: *)

(* photonNI.wl -- potentials, tortoise coordinate, Riccati-WKB boundary data and an
   arbitrary-precision integrator for the homogeneous Zerilli / Regge-Wheeler equations.

   Paper:        Sec. III D ("Numerical check") of "Gravitational radiation from a photon on
                 the light ring" (E. Barausse).

   Conventions:  M = 1, fluxes per E^2, single m > 0, j = l - m (see photonLR.wl).
                 The master equation is d^2 Psi/drstar^2 + (omega^2 - V) Psi = 0 with
                 rstar = r + 2 Log[r/2 - 1].

   Defines:
     f[r], rstar[r]                   1 - 2M/r and the tortoise coordinate;
     VRW, VZ, Vpot[l, parity, r]      Regge-Wheeler and Zerilli potentials;
     QWKB[l, parity, w, sign]         local wavenumber Q(r) of the Riccati-iterated WKB solution
                                      Psi = Exp[I Int Q drstar] (sign = +1 outgoing, -1 ingoing);
                                      three iterations of Q -> sign Sqrt[p^2 + I f Q'] starting
                                      from p = Sqrt[w^2 - V]. Used to set the boundary data of the
                                      purely ingoing solution at rstar = -40 and of the purely
                                      outgoing one at r = 40 (the normalization of each solution
                                      is fixed to unit amplitude at its own boundary);
     photonRadialNI[l, m, wp]         Psi_in, Psi_up and their rstar derivatives at r0 = 3, and the
                                      Wronskian, with WorkingPrecision wp (explicit Runge-Kutta
                                      NDSolve at arbitrary precision; slow, used for validation);
     photonModeNI[l, m, wp]           the corresponding fluxes <|"I", "H"|>.

   Loads:        photonLR.wl (source and flux formulae), from the directory of this file.
   Usage:        a package, not a script. The production run (l up to 12800,
                 photon_big_results.m) uses the machine-precision variant photonModeMP of
                 photonMP.wl, which relies on the definitions of this file. *)

schwDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
Get[FileNameJoin[{schwDir, "photonLR.wl"}]];

(* metric function and potentials (M = 1) *)
f[r_] := 1 - 2/r;
VRW[l_, r_] := f[r] (l (l + 1)/r^2 - 6/r^3);
VZ[l_, r_] := With[{lam = (l - 1) (l + 2)/2},
  f[r] (2 lam^2 (lam + 1) r^3 + 6 lam^2 r^2 + 18 lam r + 18)/(r^3 (lam r + 3)^2)];
Vpot[l_, parity_, r_] := If[parity === "Zerilli", VZ[l, r], VRW[l, r]];

(* Riccati-iterated WKB phase function Q(r): Psi = Exp[I Int Q drstar], sign = +1 outgoing, -1 ingoing *)
QWKB[l_, parity_, w_, sign_, niter_: 3] := Module[{r, p, Q},
  p = Sqrt[w^2 - Vpot[l, parity, r]];
  Q = sign p;
  Do[Q = sign Sqrt[p^2 + I f[r] D[Q, r]], {niter}];
  Function @@ {r, Q}];

rstar[r_] := r + 2 Log[r/2 - 1];

(* Arbitrary-precision integration (WorkingPrecision wp) of the two homogeneous solutions:
   Psi_in from rstar = rsA (purely ingoing WKB data) up to r0 = 3, Psi_up from r = rB (purely
   outgoing WKB data) down to r0 = 3. Returns the values and rstar derivatives at r0 and the
   Wronskian W = Psi_in Psi_up' - Psi_up Psi_in' (which equals 2 I omega A_inc, A_inc the
   incident amplitude of Psi_in at infinity, in the unit-transmission normalization). *)
photonRadialNI[l_, m_, wp_: 30, rsA_: -40, rB_: 40] := Module[{parity, w, Qin, Qup, rA, psiA, dpsiA,
    psiB, dpsiB, solIn, solUp, V, eqs, ampA, ampB, rsB, r, Psi, dPsi, rr, PsiIn, dPsiIn, PsiUp, dPsiUp, W},
  parity = If[EvenQ[l + m], "Zerilli", "ReggeWheeler"];
  w = omegaLR[m, wp + 10];
  V = Vpot[l, parity, #] &;

  (* horizon-side starting radius rA from rsA: solve rstar[r] == rsA *)
  rA = r /. FindRoot[rstar[r] == rsA, {r, 2 + 2 Exp[(rsA - 2)/2], 2 + 10^-30, 3},
    WorkingPrecision -> wp + 10, AccuracyGoal -> wp + 5];
  Qin = QWKB[l, parity, w, -1];
  Qup = QWKB[l, parity, w, +1];

  (* amplitude normalizations: |Psi_in| -> 1 at the horizon, |Psi_up| -> 1 at infinity *)
  ampA = Exp[-NIntegrate[Im[Qin[rr]]/f[rr], {rr, 2, rA},
    WorkingPrecision -> wp + 5, PrecisionGoal -> wp - 5, MaxRecursion -> 30]];
  ampB = Exp[NIntegrate[Im[Qup[rr]]/f[rr], {rr, rB, Infinity},
    WorkingPrecision -> wp + 5, PrecisionGoal -> wp - 5, MaxRecursion -> 30]];
  psiA = ampA; dpsiA = I Qin[rA] ampA;  (* phase set to zero at rA, irrelevant *)
  psiB = ampB; dpsiB = I Qup[rB] ampB;
  rsB = rstar[rB];

  (* master equation in rstar, with r(rstar) integrated alongside *)
  eqs = {Psi''[rr] + (w^2 - V[r[rr]]) Psi[rr] == 0, r'[rr] == f[r[rr]]};
  solIn = NDSolve[Join[eqs, {Psi[rsA] == psiA, Psi'[rsA] == dpsiA, r[rsA] == rA}], {Psi, r}, {rr, rsA, rstar[3]},
     WorkingPrecision -> wp, PrecisionGoal -> wp - 8, AccuracyGoal -> wp - 8, MaxSteps -> 10^7,
     Method -> "ExplicitRungeKutta", InterpolationOrder -> All][[1]];
  solUp = NDSolve[Join[eqs, {Psi[rsB] == psiB, Psi'[rsB] == dpsiB, r[rsB] == rB}], {Psi, r}, {rr, rstar[3], rsB},
     WorkingPrecision -> wp, PrecisionGoal -> wp - 8, AccuracyGoal -> wp - 8, MaxSteps -> 10^7,
     Method -> "ExplicitRungeKutta", InterpolationOrder -> All][[1]];

  (* values at r0 = 3 and Wronskian *)
  PsiIn = Psi[rstar[3]] /. solIn; dPsiIn = Psi'[rstar[3]] /. solIn;
  PsiUp = Psi[rstar[3]] /. solUp; dPsiUp = Psi'[rstar[3]] /. solUp;
  W = PsiIn dPsiUp - PsiUp dPsiIn;
  <|"PsiIn" -> PsiIn, "dPsiIn" -> dPsiIn, "PsiUp" -> PsiUp, "dPsiUp" -> dPsiUp, "W" -> W, "omega" -> w,
    "parity" -> parity|>];

(* fluxes from photonRadialNI; the jumps are converted to the tortoise coordinate,
   [dPsi/drstar] = f(r0) [dPsi/dr] *)
photonModeNI[l_, m_, wp_: 30, rsA_: -40, rB_: 40] := Module[{R = photonRadialNI[l, m, wp, rsA, rB],
    J = photonJumps[l, m], P, dP, ZI, ZH, w},
  P = J["deltaPsi"]; dP = f[3] J["deltadPsidr"];
  ZI = (R["PsiIn"] dP - P R["dPsiIn"])/R["W"];
  ZH = (R["PsiUp"] dP - P R["dPsiUp"])/R["W"];
  <|"l" -> l, "m" -> m, "omega" -> R["omega"], "Z" -> <|"I" -> ZI, "H" -> ZH|>,
    "Flux" -> photonFlux[l, m, <|"I" -> ZI, "H" -> ZH|>],
    "AbsAinc2" -> Abs[R["W"]/(2 I R["omega"])]^2, "R" -> R|>];
