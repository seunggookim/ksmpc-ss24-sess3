function Info = niftiinfogz(filename)
%NIFTIINFOGZ Unzip NII.GZ files to temp directories (works with par-for)
%
% Info = niftiinfogz(Filename)
%
% (cc) 2021-2023, sgKIM.

narginchk(1,1);
Filename = matlab.images.internal.stringToChar(filename);

if strcmp(Filename(end-2:end), '.gz')
  [~, Fname, ~] = fileparts(Filename);
  TempDir = tempname;
  gunzip(Filename, TempDir);
  TempFile = fullfile(TempDir, Fname);
  c = onCleanup(@() rmdir(TempDir, 's')); % the temp dir will be deleted whether the function exits normally or terminated by an error call.
  Info = niftiinfo(TempFile);
else
  Info = niftiinfo(filename);
end
end
