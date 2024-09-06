function h = plotridgelines(CellArray, ColorMap, vHopScale, Bndry, Scale, Markers)
% plotridgelines(CellArray, ColorMap, vHopScale, Bndry, Scale)

if not(exist('ColorMap','var')) || isempty(ColorMap)
  ColorMap = get_colormap(numel(CellArray),1);
end

if not(exist('vHopScale','var')) || isempty(vHopScale)
  vHopScale = 1;
end

if ~exist('Scale','var') || isempty(Bndry)
  Scale = ones(numel(CellArray), 1);
end

h = [];
hold on
for k = numel(CellArray):-1:1
  if exist('Bndry','var') && not(isempty(Bndry))
    [fout0,xout] = ksdensity(CellArray{k}, BoundaryCorrection='reflection', Function='pdf', Support=Bndry);
  else
    [fout0,xout] = ksdensity(CellArray{k}, BoundaryCorrection='reflection', Function='pdf');
  end

  if numel(CellArray{k})>1
    fout0 = fout0./max(fout0)*Scale(k);
    plot(xout, fout0+k*vHopScale, LineWidth=2, Color='k')
    h = [patch(xout, fout0+k*vHopScale, ColorMap(k,:), LineStyle='none'), h];
    [~,~,ci,~] = ttest(CellArray{k});
    idx = find(xout >= ci(1) & xout <= ci(2));
    if not(isempty(idx))
      patch([xout(idx(1)), xout(idx), xout(idx(end))], ...
        [0 fout0(idx) 0]+k*vHopScale, 'w', FaceAlpha=.75, LineStyle='none')
    end   
    if exist('Markers','var') && not(isempty(Markers))
      idx = dsearchn(reshape(xout,[],1), Markers(k));
      scatter(Markers(k), fout0(idx)+k*vHopScale, 80, 'o', ...
        MarkerFaceColor='w', MarkerEdgeColor='k', LineWidth=2);
    end
    
  end
end
hold off

if not(nargout)
  clear h
end

end
