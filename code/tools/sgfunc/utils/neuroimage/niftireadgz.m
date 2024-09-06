function [V,Info] = niftireadgz(varargin)
%NIFTIREADGZ Unzip NII.GZ files to temp directories (works with par-for)
%
% [V, Info] = NIFTIREADGZ(Filename)
% [V, Info] = NIFTIREADGZ(Info)
%
% (cc) 2021-2023, sgKIM.

narginchk(1,2);
varargin = matlab.images.internal.stringToChar(varargin);
if isstruct(varargin{1}) % the info syntax.
  infoIn = varargin{1};
  if isfield(infoIn, 'Filename')
    Filename = infoIn.Filename;
  else
    error(message('images:nifti:invalidStructSpecified'))
  end
else
  Filename = varargin{1};
end

if strcmp(Filename(end-2:end), '.gz')
  [~, Fname, ~] = fileparts(Filename);
  TempDir = tempname;
  gunzip(Filename, TempDir);
  TempFile = fullfile(TempDir, Fname);
  c = onCleanup(@() rmdir(TempDir, 's')); % the temp dir will be deleted whether the function exits normally or terminated by an error call.
  V = niftiread(TempFile);
  Info = niftiinfo(TempFile);
else
  V = niftiread(Filename);
  Info = niftiinfo(Filename);
end

end
