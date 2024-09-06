function [F,ker] = gaussblur(f,FWHM)
%function F = gaussblur(f,FWHM)
%
% f : input 1D series, 2D image or 3D volume
% FWHM: filtersize in sampling point (pixel/voxel).
%
% (c) Moo K. Chung, July 2003. March 2004.
% (cc) sgKIM, 2011.
% + now works for 1D, 2D and 3D.
% + kernel size extended to ensure the smoothness as intended by FWHM.
% + relationship between t and FWHM is corrected now.
%
% See Chung et al., 2001. A Unified Statistical Appraoch to 
% Deformation-based morphometry, NeuroImage for implementation detail.

if nargin == 0; help(mfilename); return; end

if ~FWHM
  F = f;
  ker = [];
  return
end

% ksize is the number of neighboring voxels we are averaging.
% it depends on FWHM.
% ksize is allowd to be an odd interger.
ksize = round(FWHM*2+1);

% see Chung et al.(2001) for the relationship between t and FWHM.
% FWHM = 4*sqrt(log(2))*sqrt(t)
%
% FWHM = sqrt(8*log(2)) * sigma (mathworld.Wolfram.com)
% sigma = FWHM / sqrt(8*log(2)) 
% t = (sigma^2)/2
t = (FWHM^2) / (16*log(2)); 

d = round(ksize/2);

% we truncate and normalize the Gaussian kernel for fast computation
% 1D, 2D, 3D version of kerenl is slightly different
if ndims(f)==1
domain_x=(1:ksize)-d;
ker=exp(-(domain_x.^2)/(4*t)) /sqrt(t/2);  % not in exponential, so doesn't matter due to normalization

elseif ndims(f)==2
domain_x=kron(ones(ksize,1),[1:ksize])-d;
domain_y=kron(ones(1,ksize),[1:ksize]')-d;
ker=exp(-(domain_x.^2 +domain_y.^2)/(4*t)) /t;

elseif ndims(f)==3
domain_x = repmat ( kron(ones(ksize,1),[1:ksize])-d, [1,1,ksize]);
domain_y = repmat ( kron(ones(1,ksize),[1:ksize]')-d, [1,1,ksize]);
domain_z = reshape( repmat ( kron(ones(ksize,1),[1:ksize])-d, [ksize, 1, 1]), [ksize,ksize,ksize]);
ker=exp(-(domain_x.^2 + domain_y.^2 +domain_z.^2)/(4*t)) /t^(3/2);

end;
ker=ker/sum(ker(:));

% smoothing is done by convolution
F=convn(f,ker,'same');



