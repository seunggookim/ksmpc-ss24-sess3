function ranovatbl = compute_peta2(ranovatbl)
% computes partial eta squared for effect size for a given RANOVA table
%
% ranovatbl = compute_peta2(ranovatbl)
%
% (cc) 2011, sgKIM.

% PARTIAL ETA^2 = (SS_effect) / (SS_effect + SS_error)


irows_int = find(contains(ranovatbl.Row, '(Intercept)'));
irows_err = find(contains(ranovatbl.Row, 'Error'));
pEta2 = [];
for iblck = 1:numel(irows_int)
  SSerr = ranovatbl.SumSq(irows_err(iblck));
  for jrow = irows_int(iblck):(irows_err(iblck)-1)
    SSeff = ranovatbl.SumSq(jrow);
    pEta2 = [pEta2; SSeff/(SSeff+SSerr)];
  end
  pEta2 = [pEta2; nan];
end
ranovatbl = addvars(ranovatbl,pEta2);
end