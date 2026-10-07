(* ::Package:: *)

(* photonMP.wl -- machine-precision integrator for the photon on the Schwarzschild light ring.
   This is the code that produced photon_big_results.m (through run_big.wls), i.e. the
   numerical fluxes of Sec. III D and Fig. 1 of "Gravitational radiation from a photon on
   the light ring" (E. Barausse). Requires the definitions of photonNI.wl (which loads
   photonLR.wl); load that file first.

   photonModeMP[l, m, rsA, rB, pg] integrates the homogeneous Zerilli (l+m even) or
   Regge-Wheeler (l+m odd) equation in the tortoise coordinate at machine precision with the
   ninth-order explicit Runge-Kutta method of NDSolve (PrecisionGoal = AccuracyGoal = pg,
   default 12), twice:
     * Psi_in from rstar = rsA (default -40), starting from the purely ingoing Riccati-WKB
       solution of photonNI.wl (unit amplitude at the horizon), up to r0 = 3;
     * Psi_up from r = rB (default 40), starting from the purely outgoing WKB solution
       (unit amplitude at infinity), down to r0 = 3.
   The two solutions, their rstar derivatives and the Wronskian W = 2 I omega A_inc at r0,
   together with the jumps of photonLR.wl (photonJumpsNum, 30-digit harmonics), give the
   amplitudes Z_I, Z_H and the fluxes. The starting radius rA on the horizon side is found
   to 40 digits so that the ingoing data are not spoiled by the exponential Delta^{-s}
   behaviour there.

   Returns an association with
     "Flux"  -> <|"I" -> Edot_I, "H" -> Edot_H|>   fluxes per E^2 (M = 1, single m),
     "AbsAinc2" -> |A_inc|^2,
     "u0", "Du0" -> Psi_in(r0)/A_inc and its rstar derivative (the u(0), u'(0) of Sec. III B),
     "v0", "Dv0" -> the same for Psi_up (the v(0), v'(0) of Sec. III C),
     "P", "dP"   -> the jumps [Psi] and [dPsi/drstar] at r0.
   (photon_big_results.m stores {l, j, Edot_I, Edot_H, u0, Du0}; the extra keys serve the
   asymmetry decomposition of asym_check.wls.)

   Cost: fractions of a second at l ~ 10, of the order of a minute per mode at l = 12800. *)

photonModeMP[l_, m_, rsA_: -40, rB_: 40, pg_: 12] := Module[{parity, w, Qin, Qup, rA, psiA, dpsiA, psiB, dpsiB, solIn, solUp, V, eqs, ampA, ampB, rsB, r, Psi, rr, PsiIn, dPsiIn, PsiUp, dPsiUp, W, J, P, dP, ZI, ZH},
  parity = If[EvenQ[l + m], "Zerilli", "ReggeWheeler"];
  w = N[m/(3 Sqrt[3])];
  V = Vpot[l, parity, #] &;
  (* horizon-side starting radius from rstar[rA] == rsA, to 40 digits *)
  rA = r /. FindRoot[rstar[r] == rsA, {r, 2 + 2 Exp[(rsA - 2)/2], 2 + 10^-30, 3}, WorkingPrecision -> 40, AccuracyGoal -> 35];
  rA = N[rA, 40];
  (* Riccati-WKB boundary data: ingoing at rA, outgoing at rB, unit amplitude at horizon / infinity *)
  Qin = QWKB[l, parity, w, -1];
  Qup = QWKB[l, parity, w, +1];
  ampA = Exp[-NIntegrate[Im[Qin[rr]]/f[rr], {rr, 2, rA}, PrecisionGoal -> 12, MaxRecursion -> 30]];
  ampB = Exp[NIntegrate[Im[Qup[rr]]/f[rr], {rr, rB, Infinity}, PrecisionGoal -> 12, MaxRecursion -> 30]];
  psiA = ampA; dpsiA = I Qin[rA] ampA;
  psiB = ampB; dpsiB = I Qup[rB] ampB;
  rsB = rstar[rB];
  (* master equation in rstar, with r(rstar) integrated alongside *)
  eqs = {Psi''[rr] + (w^2 - V[r[rr]]) Psi[rr] == 0, r'[rr] == f[r[rr]]};
  solIn = NDSolve[Join[eqs, {Psi[rsA] == psiA, Psi'[rsA] == dpsiA, r[rsA] == rA}], {Psi, r}, {rr, rsA, rstar[3]},
     PrecisionGoal -> pg, AccuracyGoal -> pg, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}, InterpolationOrder -> All][[1]];
  solUp = NDSolve[Join[eqs, {Psi[rsB] == psiB, Psi'[rsB] == dpsiB, r[rsB] == rB}], {Psi, r}, {rr, rstar[3], rsB},
     PrecisionGoal -> pg, AccuracyGoal -> pg, MaxSteps -> 10^7, Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 9}, InterpolationOrder -> All][[1]];
  PsiIn = Psi[rstar[3]] /. solIn; dPsiIn = Psi'[rstar[3]] /. solIn;
  PsiUp = Psi[rstar[3]] /. solUp; dPsiUp = Psi'[rstar[3]] /. solUp;
  W = PsiIn dPsiUp - PsiUp dPsiIn;
  (* jumps at r0 = 3, converted to the tortoise coordinate; amplitudes of the inhomogeneous solution *)
  J = photonJumpsNum[l, m, 30]; P = N[J["deltaPsi"]]; dP = N[f[3] J["deltadPsidr"]];
  ZI = (PsiIn dP - P dPsiIn)/W; ZH = (PsiUp dP - P dPsiUp)/W;
  <|"l" -> l, "m" -> m, "Flux" -> photonFlux[l, m, <|"I" -> ZI, "H" -> ZH|>], "AbsAinc2" -> Abs[W/(2 I w)]^2,
    "u0" -> PsiIn/(W/(2 I w)), "Du0" -> dPsiIn/(W/(2 I w)), "v0" -> PsiUp/(W/(2 I w)), "Dv0" -> dPsiUp/(W/(2 I w)), "P" -> P, "dP" -> dP|>];
