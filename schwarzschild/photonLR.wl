(* ::Package:: *)

(* photonLR.wl -- source, jumps and fluxes for a photon on the Schwarzschild light ring.

   Paper:    Sec. III A ("The source") and Sec. III D ("Numerical check") of
             "Gravitational radiation from a photon on the light ring" (E. Barausse).

   Problem:  a massless particle on the circular null geodesic r0 = 3M of Schwarzschild, with
             energy E, angular momentum L = b E, impact parameter b = 3 Sqrt[3] M and orbital
             frequency Omega = 1/b. Each (l, m) mode of the Zerilli (even, l+m even) or
             Regge-Wheeler (odd, l+m odd) master function Psi has frequency omega = m Omega and
             is sourced by a delta function and its derivative at r0, i.e. by a jump [Psi] and a
             jump [dPsi/dr] of the homogeneous solutions across r0.

   Conventions (used by every file of this directory):
     * units M = 1; fluxes are given per E^2 (we set E = 1), i.e. Edot M^2/E^2;
     * a single m > 0 per mode; j = l - m labels the modes near the equatorial
       dominant one (the fluxes of -m are equal, so the physical flux per l is twice
       the sum over j);
     * stress-energy coefficients (photonAK), jumps (photonJumps), Wronskian
       normalization and flux formulae (photonFlux) follow the ReggeWheeler package of
       the Black Hole Perturbation Toolkit (BHPToolkit), whose circular-orbit source is
       linear in E at fixed b = L/E, so that the photon source is the E -> 0 limit at
       fixed b (equivalently E = 1, L = b);
     * omega = m/(3 Sqrt[3]); Psi' denotes d/dr unless stated otherwise.

   Defines:
     photonAK[l, m]         Zerilli A/C/F coefficients of the point-particle source at r0;
     photonJumps[l, m]      exact jumps <|"deltaPsi", "deltadPsidr"|> at r0 (d/dr);
     photonJumpsNum[l, m]   the same with the harmonics evaluated through LogGamma (needed
                            for l of order 10^3-10^4, where SphericalHarmonicY is slow or
                            loses precision);
     photonFlux[l, m, Z]    fluxes <|"I", "H"|> to infinity / into the horizon from the
                            amplitudes Z of the inhomogeneous solution;
     photonMode[l, m, opts] fluxes from the Toolkit's own homogeneous solutions
                            (ReggeWheelerRadial, MST or numerical integration): the
                            cross-check at l <= 5 quoted in Sec. III D, run and stored by
                            ../toolkit_checks/schwarzschild_toolkit.wls (requires the
                            Toolkit; the integrators of photonNI.wl / photonMP.wl do not).

   Loads:    the ReggeWheeler package of the Toolkit, because photonMode uses it.
   Usage:    a package, not a script; it is loaded by photonNI.wl. *)

<< ReggeWheeler`

r0LR = 3;             (* light-ring radius, M = 1 *)
bLR = 3 Sqrt[3];      (* impact parameter b = L/E of the circular null geodesic *)

(* exact spherical harmonics at the equator, Y_{lm}(Pi/2, 0) and its theta derivative *)
YLR[l_, m_] := SphericalHarmonicY[l, m, Pi/2, 0];
dYLR[l_, m_] := Derivative[0, 0, 1, 0][SphericalHarmonicY][l, m, Pi/2, 0];

(* Zerilli A, C, F stress-energy coefficients (Toolkit notation) for E0 = 1, L0 = b, r0 = 3;
   the phi derivatives of the harmonic give d^2 Y/dphi^2 = -m^2 Y *)
photonAK[l_, m_] := Module[{r0 = r0LR, E0 = 1, L0 = bLR, Y, dY},
  Y = Conjugate[YLR[l, m]]; dY = Conjugate[dYLR[l, m]];
  <|"EA" -> -16 Pi (r0 - 2) E0/r0^3 Y,
    "EC" -> -16 Pi (r0 - 2) L0/r0^4/(l (l + 1)) dY,
    "EF" -> -16 Pi (r0 - 2) L0^2/E0/r0^5/((l - 1) l (l + 1) (l + 2)) (l (l + 1) - 2 m^2) Y|>];

(* jumps [Psi] and [dPsi/dr] of the master function at r0 (Toolkit formulae, Eq. (5) of
   the paper); even parity (l + m even): Zerilli; odd parity: Regge-Wheeler *)
photonJumps[l_, m_] := Module[{r0 = r0LR, ak = photonAK[l, m], EA, EC, EF, n, rm2M, np6M,
    term1, dterm1, coeff, term2},
  If[EvenQ[l + m],
    EA = ak["EA"]; EF = ak["EF"];
    n = (l - 1) (l + 2); rm2M = r0 - 2; np6M = n r0 + 6;
    term1 = -r0^4 EA/2/np6M/rm2M;
    dterm1 = -r0^4 EA/2/np6M/rm2M (4/r0 - 1/rm2M - n/np6M);
    coeff = 2/r0/rm2M;
    term2 = r0^3/4 (r0^2 n (n - 2) + r0 (14 n - 36) + 96) EA/(rm2M np6M)^2 + (n + 2) r0^2/4 EF/rm2M;
    <|"deltaPsi" -> term1, "deltadPsidr" -> -coeff term1 - dterm1 + term2|>,
    EC = ak["EC"];
    <|"deltaPsi" -> r0^3/(r0 - 2) EC,
      "deltadPsidr" -> -2 r0^2/(r0 - 2)^2 EC + r0^2/(r0 - 2) EC - 3 r0^2/(r0 - 2) EC + r0^3/(r0 - 2)^2 EC|>]];

(* mode frequency omega = m Omega = m/(3 Sqrt[3]) *)
omegaLR[m_, prec_: 40] := N[m/(3 Sqrt[3]), prec];

(* amplitudes Z of the inhomogeneous solution from the Toolkit homogeneous solutions
   R = ReggeWheelerRadial[...] (R["In"] ingoing at the horizon, R["Up"] outgoing at infinity):
   Z_I = (Psi_in [dPsi] - [Psi] Psi_in')/W, Z_H the same with Psi_up, W the Wronskian *)
photonAmplitudes[l_, m_, R_] := Module[{r0 = r0LR, J = photonJumps[l, m], PsiIn, dPsiIn, PsiUp, dPsiUp,
    W, dP, P},
  PsiIn = R["In"][r0]; dPsiIn = R["In"]'[r0]; PsiUp = R["Up"][r0]; dPsiUp = R["Up"]'[r0];
  W = PsiIn dPsiUp - PsiUp dPsiIn;
  P = J["deltaPsi"]; dP = J["deltadPsidr"];
  <|"I" -> (PsiIn dP - P dPsiIn)/W, "H" -> (PsiUp dP - P dPsiUp)/W|>];

(* energy fluxes per E^2 (Toolkit normalization of the master functions) *)
photonFlux[l_, m_, Z_] := Module[{w = omegaLR[m]},
  If[EvenQ[l + m],
    <|"I" -> (l - 1) (l + 2)/(l (l + 1)) Abs[w Z["I"]]^2/(4 Pi),
      "H" -> (l - 1) (l + 2)/(l (l + 1)) Abs[w Z["H"]]^2/(4 Pi)|>,
    <|"I" -> l (l + 1)/((l - 1) (l + 2)) Abs[w Z["I"]]^2/(16 Pi),
      "H" -> l (l + 1)/((l - 1) (l + 2)) Abs[w Z["H"]]^2/(16 Pi)|>]];

(* fluxes from the Toolkit's homogeneous solutions; opts are passed to ReggeWheelerRadial
   (e.g. WorkingPrecision -> 32, Method -> "NumericalIntegration") *)
photonMode[l_, m_, opts___] := Module[{w = omegaLR[m], R, Z},
  R = ReggeWheelerRadial[2, l, w, "Potential" -> If[EvenQ[l + m], "Zerilli", "ReggeWheeler"], opts];
  Z = photonAmplitudes[l, m, R];
  <|"l" -> l, "m" -> m, "omega" -> w, "Z" -> Z, "Flux" -> photonFlux[l, m, Z]|>];

(* Harmonics at (Pi/2, 0) for m = l - j through LogGamma: numerically stable at large l.
   Y_{l,l-j}(Pi/2,0) vanishes for odd j and its theta derivative for even j; overall signs are
   irrelevant for the fluxes. *)
YLRnum[l_, m_, prec_: 30] := Module[{j = l - m},
  If[EvenQ[j],
    N[(-1)^((2 l - j)/2) Sqrt[(2 l + 1)/(4 Pi)]
      Exp[(LogGamma[j + 1] + LogGamma[2 l - j + 1])/2 - l Log[2] - LogGamma[j/2 + 1] - LogGamma[(2 l - j)/2 + 1]],
      prec], 0]];
dYLRnum[l_, m_, prec_: 30] := Module[{j = l - m},
  If[OddQ[j], N[Sqrt[j (2 l - j + 1)], prec] YLRnum[l, m + 1, prec], 0]];

(* the A/C/F coefficients with the LogGamma harmonics *)
photonAKnum[l_, m_, prec_: 30] := Module[{r0 = r0LR, E0 = 1, L0 = N[bLR, prec], Y, dY},
  Y = YLRnum[l, m, prec]; dY = dYLRnum[l, m, prec];
  <|"EA" -> -16 Pi (r0 - 2) E0/r0^3 Y,
    "EC" -> -16 Pi (r0 - 2) L0/r0^4/(l (l + 1)) dY,
    "EF" -> -16 Pi (r0 - 2) L0^2/E0/r0^5/((l - 1) l (l + 1) (l + 2)) (l (l + 1) - 2 m^2) Y|>];

(* photonJumps with the LogGamma harmonics (used by the integrators at every l) *)
photonJumpsNum[l_, m_, prec_: 30] := Block[{photonAK = photonAKnum[#1, #2, prec] &}, photonJumps[l, m]];
