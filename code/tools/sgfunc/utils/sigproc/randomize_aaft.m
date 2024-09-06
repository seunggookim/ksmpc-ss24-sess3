function [randX, rnd_theta] = randomize_aaft(X, rnd_theta)
% [randX, rnd_theta] = randomize_aaft(X, rnd_theta)
%
% X must be a column vector
%
% (cc) 2022, dr.seunggoo.kim@gmail.com
%
% REF: Theiler et al., 1992, Testing for nonlinearity in time series: the
% method of surrogate data. Physica D 58. 77-94. doi:
% https://doi.org/10.1016/0167-2789(92)90102-S

error('Implementation error: why so large differences in magnitudes?')

% 
% assert(isvector(X)&&iscolumn(X),'vector X must be a column vector!')
% 
% % Add a duplicate of the last sample for odd-sample timeseries
% P = 1;
% [N0,~] = size(X);
% if mod(N0,2)
%   X = [X; X(end,:)];
% end
% [N,~] = size(X);
% 
% % AAFT algorithm (Theiler et al., 1992):
% % Step 1: make a Gaussian time series y[t]
% Y = sort(randn([N,P]));
% 
% % Step 2: reorder Y as X (but what about the tied ranks?)
% [~,Idx] = sort(X);
% [~,RankX] = sort(Idx);
% Y = Y(RankX);
% 
% % Step 3: randomize phases of Y -> Y'
% Yprime = randomize_phase(Y);
% 
% % Step 4: reorder X as Y'
% [~,Idx] = sort(Yprime);
% [~,RankYp] = sort(Idx);
% X = sort(X);
% randX = X(RankYp);
% 
% if mod(N0,2)
%   % drop the padded zero
%   randX(end,:) = [];
% end
% 
% end
% 
% 
% %%
% function TEST()
% X = smoothdata(randn(1000,1),'Gaussian',100);
% FT = randomize_phase(X);
% AAFT = randomize_aaft(X);
% [~,P]=kstest2(X,FT) % P<0.02
% [~,P]=kstest2(X,AAFT) % P=1
% end
