(* run_timelike.wls -- fluxes from a massive particle on circular orbits close to the light ring.

   Paper:     Sec. VI A of "Gravitational radiation from a photon on the light ring" (E. Barausse):
              the numerical fluxes against which Eq. (24) is compared. The mean of the two
              fluxes agrees with it to better than 1% for gamma = 20 and 120 <= l <= 4800, to
              0.4-4% for gamma = 10 and 30 <= l <= 1200, and to better than 4% for gamma = 5 and
              25 <= l <= 125.

   Computes:  orbits r0 = 3 (1 + delta) with delta = 1/(9 gamma^2) for gamma = 5, 10, 20, where
              gamma is the leading-order value (9 delta)^(-1/2) used in the paper; the exact
              energy per unit mass of these orbits, E0 = gamma_exact with
              gamma_exact^2 = (1+3 delta)^2/(9 delta (1+delta)), is 5.055, 10.028 and 20.014
              (printed at the start of each orbit). Eq. (24) is evaluated with the same delta.
              Modes: the dominant even mode l = m, at l = {0.1, 0.3, 1, 2, 3, 5, 8, 12}/(9 delta)
              rounded (l = 2 ... 4800), which spans the cut-off l ~ 6 gamma^2.
              Integrator: circModeMP of timelikeMP.wl (same ninth-order Runge-Kutta scheme as the
              photon run).
              The printed "pred(l)" is the leading-order prediction of Eq. (24),
              (kappa_0/l) (eta + 1/2)^2 F(eta)/F(1/2) with eta = 1/2 + 3 l delta/2.
   Inputs:    timelikeMP.wl, from the directory of this script.
   Outputs:   ../timelike_results.m, rows {r0, l, Edot_I/E^2, Edot_H/E^2} (M = 1, single m).
   Run:       from this directory: wolframscript -file run_timelike.wls (or WolframKernel -noprompt <).
   Cost:      of the order of an hour (the largest l are 1200-4800). *)

schwDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
Get[FileNameJoin[{schwDir, "timelikeMP.wl"}]];
outFile = FileNameJoin[{schwDir, "..", "timelike_results.m"}];

(* analytic prediction, m = l even sector, leading order:
   Edot_l/E^2 = (kappa0/l) (eta + 1/2)^2 F(eta)/F(1/2), eta = 1/2 + 3 l delta/2 *)
F[a_] := Exp[-Pi a/2]/(Cosh[Pi a] Abs[Gamma[3/4 + I a/2]]^2);
kappa0 = Sqrt[Pi] Sech[Pi/2]/(9 Exp[Pi/4] Abs[Gamma[3/4 + I/4]]^2);
pred[l_, d_] := Module[{a = 1/2 + 3 l d/2}, (kappa0/l) (a + 1/2)^2 F[a]/F[1/2]];

res = {};
Do[
  d = (r0 - 3)/3; gam2 = (1 + 3 d)^2/(9 d (1 + d));       (* delta and the exact gamma^2 = E0^2 *)
  Print["r0=", r0, " delta=", d, " gamma=", Sqrt[gam2]];
  Do[
    md = circModeMP[l, l, r0];
    e2 = md["E0"]^2;                                        (* fluxes per mu^2 -> per E^2 *)
    AppendTo[res, {r0, l, md["Flux"]["I"]/e2, md["Flux"]["H"]/e2}];
    Print["  l=", l, " l*Edot_I/E^2=", l md["Flux"]["I"]/e2, " l*Edot_H/E^2=", l md["Flux"]["H"]/e2,
      "  mean=", l (md["Flux"]["I"] + md["Flux"]["H"])/(2 e2), "  pred(l)=", l pred[l, d],
      "  exp(-pi l d 3/2)=", Exp[-3 Pi l d/2]],
    {l, Round[{0.1, 0.3, 1, 2, 3, 5, 8, 12}/(9 d)]}];
  , {r0, {3 (1 + 1/(9 25.)), 3 (1 + 1/(9 100.)), 3 (1 + 1/(9 400.))}}];   (* gamma = 5, 10, 20 *)
Put[res, outFile];
Print["DONE"];
