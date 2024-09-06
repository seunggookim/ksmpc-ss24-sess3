function niftiwritegz(V, filename, varargin)
% NIFTIWRITEGZ MATLAB's gzip cannot handle '~', so expand it then run it.
% IGNORE when filename includes extension '.nii' or '.nii.gz'
%
% [] = NIFTIWRITEGZ(Img, Filename)
% [] = NIFTIWRITEGZ(Img, Filename, Info)
% [] = NIFTIWRITEGZ(Img, Filename, Info, ...)
%
% (cc) 2021-2023, sgKIM.

[p1, f1, e1] = myfileparts(filename);
Idx = find(cellfun(@(x) strcmpi(x, 'Compressed'), varargin));
if ~isempty(e1)     % if extension is defind:
  if ~isempty(Idx)  % and if COMPRESSED is defined:
    if strcmpi('.nii.gz', e1) ~= varargin{Idx+1}
      warning([...
        'COMPRESSED=%i is not consistent with the extension "%s": ',...
        'COMPRESSED is corrected according to the extension'], ...
        varargin{Idx+1}, e1)
    end
    varargin{Idx+1} = strcmpi('.nii.gz', e1); % make it consistent
  else              % if COMPRESSED is not defined:
    varargin = [varargin, 'Compressed', strcmpi('.nii.gz', e1)];
  end
end

niftiwrite(V, fullfile(p1, f1), varargin{:})
end
