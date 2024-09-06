function myeeg_plotchan(Colormap, Chanlocs)

topoplot([], Chanlocs, 'style','blank', 'whitebk','on');
h = get(gca, 'Children');
h(1).Visible = 'off';
patch(h(5).XData, h(5).YData, 0, 'FaceColor','w');
hold on; scatter(h(1).XData, h(1).YData, h(1).MarkerSize*3, Colormap, 'filled'); hold off

end