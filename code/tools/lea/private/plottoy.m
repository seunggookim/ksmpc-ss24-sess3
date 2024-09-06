function plottoy(X, Y, Mdl, Rnd, Job)
% plottoy(X, Y, Mdl, Rnd, Job)

Data = conformdata(X, Y, Job);

set(gcf, Colormap = flipud(brewermap(256,'Spectral')) );
AxesPos = axeslayout([3 3]);

axespos(AxesPos,1); imagesc(X{1}); colorbar; 
title(sprintf('X_1 | TempGaussWin=%i smp', Job.TempGaussWin))
cL = clim(); clim([-max(abs(cL)), +max(abs(cL))])
xlabel('Features#'); ylabel('Sample#');

axespos(AxesPos,2); imagesc(Y{1}); colorbar;
cL = clim(); clim([-max(abs(cL)), +max(abs(cL))])
title(sprintf('Y_1 | EffectSize=%g',Job.EffectSize))
xlabel('Response#'); ylabel('Sample#')

axespos(AxesPos,4); imagesc(Data(1).X); colorbar; title('Xc_1')
cL = clim(); clim([-max(abs(cL)), +max(abs(cL))])
xlabel('Predictor#'); ylabel('Sample#');

axespos(AxesPos,5); imagesc(Data(1).Y); colorbar;
cL = clim(); clim([-max(abs(cL)), +max(abs(cL))])
title(sprintf('Yc_1 | EffectSize=%g',Job.EffectSize))
xlabel('Response#'); ylabel('Sample#')

axespos(AxesPos,3); imagesc(squeeze(mean(Mdl.Bhat,1)));
cL = clim(); clim([-max(abs(cL)), +max(abs(cL))])
colorbar; title('mean B-hat [w/ intercept]')
xlabel('Response#'); ylabel('Predictor#')

axespos(AxesPos,6); imagesc(log10(Mdl.Lopt)); colorbar; title('log_{10}\lambda')
cL = clim(); clim([-max(abs(cL)), +max(abs(cL))])
xlabel('Response#'); ylabel('Fold#')

axespos(AxesPos,7); imagesc(Mdl.Acc); colorbar; title('Acc [r]')
clim([-1 1]); xlabel('Response#'); ylabel('Fold#')

axespos(AxesPos,8);
nResp = size(Rnd.AccRnd,2);
X_ = cellfun(@(x) Rnd.AccRnd(:,x), num2cell(1:nResp), 'Uni',0);
plotridgelines(X_, brewermap(nResp,'Set2'), 1, [], [], mean(Mdl.Acc,1))
ylabel('Response#'); xlabel('Pred. Acc. [r]'); xlim([-1 1])
title(sprintf('Obs. vs. null pred. acc. | k=%i', Job.nRands))

axespos(AxesPos,9);
bar(-log10(Rnd.PvalFdr)); xlabel('Response#'); ylabel('-log_{10}(P_{fdr})')
yline(-log10(0.01), color=[0.8500    0.3250    0.0980], LineWidth=2)
title('FDR-adjusted P-values')
end
