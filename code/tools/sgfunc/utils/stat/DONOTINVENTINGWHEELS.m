%% Let's use STAT & ML TOOLBOX (introduced in R2012a)
clear;clc
rng(sum(double('Seung-Goo Kim')))
X = rescalecol(randn(100,10));
% X = randn(100,10);
% X(:,2) = X(:,2) + X(:,1)*0.5;
beta = randn(10,1)

rng('shuffle')
Y = X*beta + randn(100,1)*1 + X(:,1).*X(:,2)*0;

mdl = fitlm(X,Y,'interactions')

tic
mdl=fitlm(X,Y);
mdl.Residuals.Raw;
toc

% tic
% [y1,~,y2]=residual(Y+1000,X);
% toc
% figure;
% subplot(121)
% histfit(y1)
% subplot(122)
% histfit(y2)


%%
clf
subplot(441)
stat_desc(Y)

subplot(445)
h1=plotEffects(mdl);
hold on
h2=scatter(range(X)'.*beta, 1:10,'rx');
legend([h1(1), h2],{'Est','True'},'box','on','location','northoutside')

subplot(446)
plotInteraction(mdl,'x1','x2')
subplot(447)
plotAdded(mdl)
title(sprintf('adj R^2 = %.2f',mdl.Rsquared.Adjusted))
subplot(448)
plotResiduals(mdl)

subplot(449)
h = plotAdded(mdl,1+[8],'markerEdgeColor','m');
% h(3).LineStyle = '-';
subplot(4,4,10)
h = plotAdded(mdl,1+[1 8],'markerEdgeColor','m');

%% Mixed effect models: repeated measures model
clear; clf; clc
load repeatedmeas
rm = fitrm(between,'y1-y8 ~ Group*Gender+Age+IQ','WithinDesign',within)
subplot(441);
imagesc(table2array(between(:,[1:8]+4)))
subplot(442)
h = plot(rm);
cmap = turbo(numel(h));
for i=1:numel(h)
  h(i).Color = cmap(i,:);
end
rtbl = ranova(rm,'WithinModel','w1+w2+w1*w2')

%% Nonlinear regression
clear; clf; clc;
load carbig
X = [Horsepower,Weight];
subplot(441)
imagesc(zscore(X))
y = MPG;
subplot(442)
histfit(y)
modelfun = @(b,x) b(1) + b(2)*x(:,1).^b(3) + b(4)*x(:,2).^b(5);
beta0 = [-50 200 -1 200 -1];
mdl = fitnlm(X, y, modelfun, beta0)
subplot(443)
plotAdded(fitlm(y,predict(mdl)))
subplot(444)
plotResiduals(mdl)

% new figure:
plotSlice(mdl)