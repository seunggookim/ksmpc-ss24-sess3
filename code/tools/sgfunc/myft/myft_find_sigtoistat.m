function [XLIMno,XLIM,PEAK]=myft_find_sigtoistat(stat)
if isfield(stat,'posclusters') || isfield(stat,'negclusters')
  % find time interval of significant clusters:
  XLIM={};PEAK=[];
  % find significant positive clusters
  if isfield(stat,'posclusters')
    clus=stat.posclusters;
    idx_clus=find([clus.prob]<stat.cfg.alpha);
    for j=1:numel(idx_clus)
      mask=stat.posclusterslabelmat==idx_clus(j);
      trms=rms(denan(stat.stat.*mask),1);
      [~,idx_peak]=max(trms);
      PEAK=[PEAK stat.time(idx_peak)];
      ind_time=any(mask,1);
      xlim=stat.time(ind_time);
      xlim=xlim([1 end]);
      XLIM=[XLIM xlim];
    end
  end
  % find significant negative clusters
  if isfield(stat,'negclusters')
    clus=stat.negclusters;
    idx_clus=find([clus.prob]<stat.cfg.alpha);
    for j=1:numel(idx_clus)
      %ind_time=any(stat.negclusterslabelmat==idx_clus(j),1);
      mask=stat.negclusterslabelmat==idx_clus(j);
      trms=rms(denan(stat.stat.*mask),1);
      [~,idx_peak]=max(trms);
      PEAK=[PEAK stat.time(idx_peak)];
      ind_time=any(mask,1);
      xlim=stat.time(ind_time);
      xlim=xlim([1 end]);
      XLIM=[XLIM xlim];
    end
  end
  [PEAK,idx_peak]=sort(PEAK); % reorder time points
  XLIM=XLIM(idx_peak);
else
  % find RMS peaks
  maskedT=stat.stat.*stat.mask;
  y=rms(denan(maskedT));
  [pks,idx_pks]=findpeaks(y);
  PEAK=stat.time(idx_pks);
end

% find overlaps and collapse if overlap is >90%
if numel(XLIM)>1
  XLIMno={};
  XLIMno{1}=XLIM{1};
  j=1;
  for i=2:numel(XLIM)
    a=XLIMno{j};
    b=XLIM{i};
    if a(1)>b(1)
      b=XLIMno{j};
      a=XLIM{i};
    end
    ovlp=(a(2)-b(1))/(b(2)-a(1))
    if ovlp>0.5
      XLIMno{j}=[a(1) max([a(2) b(2)])];
    else
      j=j+1;
      XLIMno{j}=XLIM{i};
    end
  end
else
  XLIMno=XLIM;
end

end