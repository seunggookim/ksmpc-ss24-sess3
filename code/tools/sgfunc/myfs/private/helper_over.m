function h = helper_over(Ha, Vi, u, w, cfg)

if isequal(cfg.thres, [0 0]) % no thresholding
  Vi(Vi==0) = nan;
  if size(Vi,3) == 3 % an RGB-slice
    h = pcolor(Ha, u.axis, w.axis, Vi(:,:,1)*0);
    h.CData = Vi;
    h.LineStyle = 'none';
  else
    h = pcolor(Ha, u.axis, w.axis, Vi);
    h.LineStyle = 'none';
  end
else
  % suprathreshold image:
  Vsupra = Vi;
  Vsupra(Vsupra==0) = nan;
  Vsupra(cfg.thres(1)<Vi & Vi<cfg.thres(2)) = nan;
  h = pcolor(Ha, u.axis, w.axis, Vsupra);
  h.LineStyle = 'none';

  % subthreshold image:
  if cfg.subthres
    hold on
    Vsub = Vi;
    Vsub(Vsub==0) = nan;
    Vsub(Vi<cfg.thres(1) | cfg.thres(2)<Vi) = nan;
    h2 = pcolor(Ha, u.axis, w.axis, Vsub);
    h2.LineStyle = 'none';
    h2.FaceAlpha = 0.5;
    hold off
  end
end
set(Ha,'color','none')
end
