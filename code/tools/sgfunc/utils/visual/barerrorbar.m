function h = barerrorbar(x, meany, erry)%, facecolor, edgecolor, width, errbar_bidir)
% h = barerrorbar(x, meany, erry)
% (cc) dr.seunggoo.kim@gmail.com

hBar = bar(x, meany);
hLine = [];
for i = 1:numel(hBar)
  for j = 1:numel(hBar(i).XEndPoints)
    x = hBar(i).XEndPoints(j);
    y = hBar(i).YEndPoints(j);
    hLine = [hLine, line([x, x], [y-erry(i,j), y+erry(i,j)], Color=hBar(i).FaceColor, LineWidth=2)];
  end
end

end