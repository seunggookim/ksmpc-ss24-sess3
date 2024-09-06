function [w,gi,lag,ri] = myspm_est_res(SPM)
% code taken from spm_est_V.m
cd(SPM.swd)
dt = SPM.xY.RT;
X  = SPM.xX.X;
X  = [X SPM.xX.K.X0];                         % add low-freq cosine bases
X  = [spm_conv(randn(size(X,1),4),8/dt,0) X]; % add smooth random reg
m  = size(X,1);
R  = speye(m,m) - X*spm_pinv(X);              % residual forming matrix

% get data from significant voxels
%--------------------------------------------------------------------------
N     = 4000;                               % number of voxels
Vspm  = SPM.xCon(1).Vspm;                   % get first SPM
XYZ   = SPM.xVol.XYZ;
F     = spm_sample_vol(Vspm,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
[F,i] = sort(F,2,'descend');
XYZ   = XYZ(:,i(1:16000));                  % voxels for t-test

% get data and covariance
%--------------------------------------------------------------------------
Y     = spm_get_data(SPM.xY.VY,XYZ);
m     = size(Y,1);                          % number of scans
% Y   = spm_null_data(Y,SPM);               % uncomment for null data
 
% Data covariance
%--------------------------------------------------------------------------
C     = cov(Y(:,1:N)');
% C  = SPM.xVi.Cy;                              % sample covariance
C  = C*trace(R*R')/trace(R*C*R');             % scaling

% % covariance components (a mixture of exponentials)
% %==========================================================================
% dt    = SPM.xY.RT;                            % TR (seconds)
% T     = (0:(m - 1))*dt;                       % time
% d     = 2.^(floor(log2(dt/4)):log2(64));      % time constants (seconds)
% QQ    = {};                                   % dictionary of components
% for i = 1:length(d)
%     for j = 0:1
%         QQ{end + 1} = toeplitz((T.^j).*exp(-T/d(i)));
%     end
% end
% Q     = QQ(1:end);                            % full (exactly)
% [V,h] = spm_reml(C,X,Q,1,1,0,4);              % ReML to est V
V  = SPM.xVi.V;                               % intrinsic cov after HPF
W  = spm_inv(spm_sqrtm(V));                   % W*W' = V^-1
W  = W*sqrt(trace(R*R')/trace(R*W*C*W'*R'));  % scaling

% Power spectra:
S  = spm_sqrtm(C);
g  = [  sum(abs(fft(full(R)).^2),2)];         % residual forming matrix
g  = [g sum(abs(fft(full(R*S)).^2),2)];       % residuals unwhitened
g  = [g sum(abs(fft(full(R*W*S)).^2),2)] ;    % residuals whitened
i  = floor(2:m/2);
w  = (1:length(i))/(2*length(i)*dt);
gi = g(i,:);

% Autocovariance:
r     = [];
for i = 1:size(g,2)
  f      = ifft(g(:,i));
  r(:,i) = real(fftshift(f));
end
lag   = -32:32;
iseven = @(x) ~mod(x,2);
i     = lag + fix(m/2) + iseven(m);
ri    = r(i,:);
lag   = lag*dt;

if ~nargout
  figure
  subplot(121)
  plot(w,gi)
  legend({'ideal','unwhitened','whitened'})
  xlabel('Freq [Hz]'); ylabel('Power')
  title('Spectral density')
  subplot(122)
  plot(lag,ri)
  xlabel('Lag [s]'); ylabel('Cov');
  title('Autocovariance func')
end

end