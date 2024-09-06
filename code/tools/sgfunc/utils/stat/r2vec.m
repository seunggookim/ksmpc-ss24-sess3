function r2 = r2vec(Y,Yhat)
% r2 = r2vec(y,yhat)
%
% computes R2 for each column vector of Y (NxP) given Yhat (NxP):
% Y = [y1, y2, ..., yp]
% y1 = [y11, y12, ..., y1n]'
% 
%             sum( (yi - yi_hat)^2 )
% R2_i = 1 - ------------------------
%             sum( (yi - yi_bar)^2 )
%
% where yi_bar is a mean of yi, yi_hat is a prediction for yi.
%
% (cc) 2019, sgKIM. solleo@gmail.com

if ~(all(size(Y)==size(Yhat)))
  error('Y and Yhat should be in the same dimension!')
end
Ybar = ones(size(Y)).*mean(Yhat,1);
r2 = 1 - sum( (Y-Yhat).^2 )./sum( (Y-Ybar).^2 );
end
