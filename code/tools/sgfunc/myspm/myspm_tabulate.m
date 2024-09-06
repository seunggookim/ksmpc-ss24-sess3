function myspm_tabulate (tmap, clustermap, fname_xls, isfslprob)

% if ~exist('mask','var')
%   mask = 0.*tmap + 1;
% end
% clustermap = clustermap .* mask;

if ~exist('isfslprob')
isfslprob=0;
end

idx    = unique(clustermap(:));
idx(idx==0) = [];

if ~numel(idx)
warning('no sig.');
return
end

if ~exist('fname_tsv','var')
fname_xls=['table_', datestr(now,'YYMMDD_hhmmss') '.xls'];
end

for j=1:numel(idx) % for each cluster
tmap1 = tmap.*(clustermap == idx(j));
% find the name of the cluster
if isfslprob
[cst, cst_all] =myfsl_atlasquery(find3(tmap1),1);
clusname=[cst.name,' (',num2str(cst.prob),'%)'];
else
cst=myspm_NMatlas(find3(tmap1),1);
clusname=[cst.name,' (',num2str(cst.prob),'%)'];
end
disp(['#',num2str(j),':clus: ',clusname]);

% find the peak
tvals = tmap(clustermap == idx(j));
if mean(tvals) >0
peak_ijk = find3 (tmap1 == max(tvals));
elseif mean(tvals) <0
peak_ijk = find3 (tmap1 == min(tvals));
end

% find the name of peak
if isfslprob
[st, st_all]=myfsl_atlasquery(peak_ijk,1);
peakname=[st.name,' (',num2str(st.prob),'%)'];
else
st=myspm_NMatlas(peak_ijk,1);
peakname=[st.name];
end
disp(['#',num2str(j),':peak: ',peakname,', T=',num2str(tmap1(peak_ijk(1), peak_ijk(2), peak_ijk(3)))]);

end


disp(['Saving results in: ',fname_xls]);


end
