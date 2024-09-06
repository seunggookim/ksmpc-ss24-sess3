function JOB = myspm_mycon (JOB)
% JOB = myspm_mycon (JOB)
% computes contrast images with a given name con_*
%
% JOB requires:
%  .dir_glm  [1xN string]
%  .conname  [1xK cell]   names of K contrasts
%  .convec   [KxP double] contrast vectors for K contrasts, when the number
%                        of betas' is larger than P, zeros will be padded.
%
% (cc) 2016, sgKIM

dir0=pwd;
cd(JOB.dir_glm)
convec = JOB.convec;
if ~isfield(JOB,'overwrite'), overwrite=0; else overwrite=JOB.overwrite; end
if exist('beta_0001.img','file') && ~exist('beta_0001.nii','file')
 ext1='.img';
elseif ~exist('beta_0001.img','file') && exist('beta_0001.nii','file')
 ext1='.nii';
else
 error(['No beta files found in ',JOB.dir_glm])
end
% set contrast names
K=size(convec,1);
if ~isfield(JOB,'conname')
 warning(['No contrast name is given. Will save as con_0001',ext1]);
 for k=1:K
  conname{k} = ['con_',pad(k,4)];
 end
end
conname = JOB.conname;
alldone=1;
for k=1:K
 alldone= alldone*exist([JOB.dir_glm,'/con_',conname{k},ext1],'file');
end
if alldone && ~overwrite, return; end

% find variable names from SPM.mat
varnames=[];
if isfield(JOB,'varnames')
 varnames=JOB.varnames;
else
 if exist('variablenames.mat','file')
  load('variablenames.mat','varnames');
 elseif exist('SPM.mat','file')
  D=dir('SPM.mat');
  if D.bytes/1024^3 > 0.1
   warning('SPM file too big (>0.1 GB), set "varnames" to be []');
  else
   load SPM.mat
   varnames = SPM.xX.name;
   save('variablenames.mat','varnames');
  end
 end
end

% find beta's
fnames_beta = mydir(['beta*',ext1]);
if isempty(fnames_beta)
 error(['No beta files found'])
end
B=numel(fnames_beta); % B= # of beta files
[K,P]= size(convec);  % P= # of elements of the contrast vector

% pad zeros if needed
if P<B
 warning('Padding zeros to the given contrast vector..');
 convec = [convec, zeros(K, (B-P))];
elseif P>B
 warning(['Given contrast vector has more elements than computed beta''s: ', ...
  'only elements from the left are used.']);
end

% read beta's
disp('> reading beta''s..');
hdr=load_untouch_header_only(fnames_beta{1});
BETA=zeros([hdr.dime.dim(2:4) B],'single');
h=figure;
for b=1:B
 nii=load_uns_nii(fnames_beta{b});
 vol=nii.img;
 vol(~~isnan(vol))=0;
 BETA(:,:,:,b) = vol;
end
for k=1:K
 contraststr='';
 con=BETA(:,:,:,1)*0;
 normcon=convec(k,:);
 % scaling
 idx1=normcon>0; idx2=normcon<0;
 if (sum(abs(normcon(idx1))) > 1) || (sum(abs(normcon(idx2))) > 1)
  normcon(idx1) = normcon(idx1)./sum(abs(normcon(idx1)));
  normcon(idx2) = normcon(idx2)./sum(abs(normcon(idx2)));
 end
 for b=1:B
  c1 = normcon(b);
  con = con + BETA(:,:,:,b).*c1;
  if ~isempty(varnames) && ~~c1
   Sign='-+';
   contraststr=[contraststr, ...
    Sign((c1>0)+1),num2str(abs(c1)),'*[',varnames{b},'] '];
  end
 end
 if ~isempty(varnames), disp(['> contrast: ',contraststr]); end
 clf; imageorth(con); title(['con_',conname{k},ext1]); drawnow;
 nii.img = con;
 nii.hdr.hist.descrip = ...
  ['myspm: dir_glm=',JOB.dir_glm,'; contrast=', conname{k}];
 fname_con=[JOB.dir_glm,'/con_',conname{k},ext1];
 disp(['> writing ',fname_con]);
 disp(' ')
 save_untouch_nii(nii, fname_con);
end
cd (dir0);
close(h);
end
