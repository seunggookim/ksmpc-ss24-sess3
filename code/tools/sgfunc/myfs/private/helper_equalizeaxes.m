function helper_equalizeaxes(baseaxes, overaxes)
xlim = [min(baseaxes.XLim(1),overaxes.XLim(1)), max(baseaxes.XLim(2),overaxes.XLim(2))];
ylim = [min(baseaxes.YLim(1),overaxes.YLim(1)), max(baseaxes.YLim(2),overaxes.YLim(2))];
baseaxes.XLim = xlim;
baseaxes.YLim = ylim;
overaxes.XLim = xlim;
overaxes.YLim = ylim;
end
