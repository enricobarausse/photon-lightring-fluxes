(* asym_check.wls -- decomposition of the difference between the two fluxes for j = 0.
   Supports the end of Sec. III D of "Gravitational radiation from a photon on the light ring"
   (E. Barausse): "at l = 6400 the cross term is 1.0389, against the analytic 1.0390, and the
   barrier asymmetry is -0.231".

   With u = Psi_in/A_inc, v = Psi_up/A_inc at r0 and the jumps P = [Psi], dP = [dPsi/drstar],
   the fluxes are proportional to I0 = |u dP - P u'|^2 and H0 = |v dP - P v'|^2, and
   Sqrt[l] (I0 - H0)/(I0 + H0) splits into
     A1 = Sqrt[l] (|u|^2 - |v|^2) dP^2/den          barrier asymmetry (cubic term of the barrier),
     A2 = -2 Sqrt[l] (Re[u u'*] - Re[v v'*]) dP P/den   cross term (analytic limit 1.0390),
     A3 = Sqrt[l] (|u'|^2 - |v'|^2) P^2/den         O(1/l) remainder,
   den = (|u|^2 + |v|^2) dP^2. The quantities come from photonModeMP of photonMP.wl.

   Run from this directory: wolframscript -file asym_check.wls. Cost: a few minutes
   (dominated by l = 6400). *)

schwDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];
Get[FileNameJoin[{schwDir, "photonNI.wl"}]];
Get[FileNameJoin[{schwDir, "photonMP.wl"}]];

Print["l, sqrt(l)*(I-H)/(I+H) [j=0], A1 (|u0|^2-|v0|^2 term), A2 (cross term; pred 1.0390), A3 (|Du0|^2[Psi]^2 terms, O(l^-1))"];
Do[md = photonModeMP[l, l];
  u0 = md["u0"]; Du0 = md["Du0"]; v0 = md["v0"]; Dv0 = md["Dv0"]; P = md["P"]; dP = md["dP"];
  I0 = Abs[u0 dP - P Du0]^2; H0 = Abs[v0 dP - P Dv0]^2;
  den = (Abs[u0]^2 + Abs[v0]^2) dP^2;
  A1 = Sqrt[l] (Abs[u0]^2 - Abs[v0]^2) dP^2/den;
  A2 = Sqrt[l] (-2 (Re[u0 Conjugate[Du0]] - Re[v0 Conjugate[Dv0]]) dP P)/den;
  A3 = Sqrt[l] ((Abs[Du0]^2 - Abs[Dv0]^2) P^2)/den;
  Print[{l, Sqrt[l] (I0 - H0)/(I0 + H0), A1, A2, A3}], {l, {100, 400, 1600, 6400}}];
