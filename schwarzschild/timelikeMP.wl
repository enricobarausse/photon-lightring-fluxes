(* ::Package:: *)

(* timelikeMP.wl -- massive particle on a timelike circular orbit of radius r0 (M = 1),
   same integrator as photonMP.wl. Supports Sec. VI A ("Timelike orbits near the light ring")
   of "Gravitational radiation from a photon on the light ring" (E. Barausse): the comparison
   of Eq. (timelike) with numerical fluxes for gamma = 5-20 (run_timelike.wls).
   Loads photonNI.wl (potentials, WKB data) and photonMP.wl.

   Conventions: Toolkit point-particle source with E0 = Sqrt[(r0-2)(r0^2+L0^2)/r0^3] and
   L0 = Sqrt[r0^2/(r0-3)] per unit mass mu, so the fluxes returned are per mu^2; the energy
   per unit mass E0 is returned too, so that Flux/E0^2 is the flux per E^2 used in the paper.
   Orbital frequency Omega = r0^(-3/2), omega = m Omega; r0 = 3 (1 + delta), gamma = E0.

   circJumps[l, m, r0]     jumps [Psi], [dPsi/dr] at r0 (Toolkit formulas, LogGamma harmonics);
   circModeMP[l, m, r0]    fluxes <|"I", "H"|> per mu^2 and "E0", with the ninth-order explicit
                           Runge-Kutta integration from rstar = -40 (ingoing WKB data) and from
                           r = 40 (outgoing WKB data) to r0, as in photonModeMP. *)

schwDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
Get[FileNameJoin[{schwDir, "photonNI.wl"}]];
Get[FileNameJoin[{schwDir, "photonMP.wl"}]];

circJumps[l_, m_, r0_, prec_: 30] := Module[{E0, L0, Y, dY, EA, EC, EF, n, rm2M, np6M, term1, dterm1, coeff, term2},
  L0 = N[Sqrt[r0^2/(r0 - 3)], prec]; E0 = N[Sqrt[(r0 - 2) (r0^2 + L0^2)/r0^3], prec];
  Y = YLRnum[l, m, prec]; dY = dYLRnum[l, m, prec];
  EA = -16 Pi (r0 - 2) E0/r0^3 Y;
  EC = -16 Pi (r0 - 2) L0/r0^4/(l (l + 1)) dY;
  EF = -16 Pi (r0 - 2) L0^2/E0/r0^5/((l - 1) l (l + 1) (l + 2)) (l (l + 1) - 2 m^2) Y;
  If[EvenQ[l + m],
    n = (l - 1) (l + 2); rm2M = r0 - 2; np6M = n r0 + 6;
    term1 = -r0^4 EA/2/np6M/rm2M;
    dterm1 = -r0^4 EA/2/np6M/rm2M (4/r0 - 1/rm2M - n/np6M);
    coeff = 2/r0/rm2M;
    term2 = r0^3/4 (r0^2 n (n - 2) + r0 (14 n - 36) + 96) EA/(rm2M np6M)^2 + (n + 2) r0^2/4 EF/rm2M;
    <|"deltaPsi" -> term1, "deltadPsidr" -> -coeff term1 - dterm1 + term2, "E0" -> E0|>,
    <|"deltaPsi" -> r0^3/(r0 - 2) EC,
      "deltadPsidr" -> -2 r0^2/(r0 - 2)^2 EC + r0^2/(r0 - 2) EC - 3 r0^2/(r0 - 2) EC + r0^3/(r0 - 2)^2 EC, "E0" -> E0|>]];

circModeMP[l_, m_, r0_, rsA_: -40, rB_: 40, pg_: 12] := Module[{parity, w, Qin, Qup, rA, psiA, dpsiA, psiB, dpsiB, solIn, solUp, V, eqs, ampA, ampB, rsB, r, Psi, rr, PsiIn, dPsiIn, PsiUp, dPsiUp, W, J, P, dP, ZI, ZH, rs0, flI, flH, ce},
  parity = If[EvenQ[l + m], "Zerilli", "ReggeWheeler"];
  w = N[m Sqrt[1/r0^3]];
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
  J = circJumps[l, m, r0]; P = N[J["deltaPsi"]]; dP = N[f[r0] J["deltadPsidr"]];
  ZI = (PsiIn dP - P dPsiIn)/W; ZH = (PsiUp dP - P dPsiUp)/W;
  (* flux normalisation of photonFlux, written out since omega is not the light-ring one *)
  ce = If[EvenQ[l + m], (l - 1) (l + 2)/(l (l + 1))/(4 Pi), l (l + 1)/((l - 1) (l + 2))/(16 Pi)];
  <|"l" -> l, "m" -> m, "r0" -> r0, "E0" -> J["E0"], "Flux" -> <|"I" -> ce Abs[w ZI]^2, "H" -> ce Abs[w ZH]^2|>|>];
