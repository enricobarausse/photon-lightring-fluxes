(* run_big.wls -- production run for the Schwarzschild light ring: fluxes up to l = 12800.

   Paper:     the data of Fig. 1 and of the numbers quoted in Sec. III D of "Gravitational
              radiation from a photon on the light ring" (E. Barausse) (three-term fit, kappa_j at
              l = 12800, m-summed asymmetry).

   Computes:  for each l in {10, 20, 50, ..., 12800} and j = l - m = 0..4 (the five dominant m > 0),
              photonModeMP of photonMP.wl (ninth-order explicit Runge-Kutta NDSolve, Riccati-WKB
              boundary data) gives the fluxes to infinity and into the horizon per E^2 (M = 1).
   Inputs:    photonNI.wl and photonMP.wl, from the directory of this script.
   Outputs:   ../photon_big_results.m, rows {l, j, Edot_I, Edot_H, u(0), u'(0)}, single m; the
              physical flux per l is twice the sum over j (the -m modes contribute the same). The
              file is rewritten after every mode, so a partial run is usable.
   Run:       from this directory (or with wolframscript -file from anywhere):
                wolframscript -file run_big.wls           or      WolframKernel -noprompt < run_big.wls
              Mathematica 15 with the ReggeWheeler package of the Black Hole Perturbation Toolkit
              installed (loaded by photonLR.wl; it is not used by the integration itself).
              Quick test: l values on the command line, e.g. wolframscript -file run_big.wls 10 20,
              write to photon_big_test.m in the current directory instead of the stored file.
   Cost:      the full run takes of the order of an hour (about a minute per mode at l = 12800). *)

schwDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
Get[FileNameJoin[{schwDir, "photonNI.wl"}]];
Get[FileNameJoin[{schwDir, "photonMP.wl"}]];

ls = {10, 20, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800};
outFile = FileNameJoin[{schwDir, "..", "photon_big_results.m"}];

(* quick test: l values from the command line, output in the current directory *)
If[Length[$ScriptCommandLine] > 1,
  ls = ToExpression /@ Rest[$ScriptCommandLine]; outFile = "photon_big_test.m"];

res = {};
Do[
  Do[
    t = AbsoluteTiming[md = photonModeMP[l, l - j];];
    AppendTo[res, {l, j, md["Flux"]["I"], md["Flux"]["H"], md["u0"], md["Du0"]}];
    Print["l=", l, " j=", j, " t=", NumberForm[t[[1]], 3], "  l*Edot_I=", l md["Flux"]["I"],
      "  l*Edot_H=", l md["Flux"]["H"], "  sum=", l (md["Flux"]["I"] + md["Flux"]["H"])];
    Put[res, outFile];                     (* rewritten after every mode *)
    , {j, 0, 4}];
  , {l, ls}];
Print["DONE"];
