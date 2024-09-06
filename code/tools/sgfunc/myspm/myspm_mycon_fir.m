function JOB = myspm_mycon_fir(JOB)
% JOB = myspm_mycon_fir(JOB)
%
% this does... what?
%
% JOB requires:
%  .dir_glm
%  .cntname_fir
%  .cntvec_fir
%  .num_delay
%  .TR_sec
%  .length_sec
%  .winsize
%  .num_cond
%
% (cc) 2016, sgKIM

num_delay  = JOB.num_delay;
length_sec = JOB.length_sec;
TR_sec  = JOB.TR_sec;
winsize = JOB.winsize;
postonset_sec=[];
cntMtx=[];
cntName={};
for k=1:size(JOB.cntvec_fir,1)
 if ~mod(winsize,2)
  winsize=winsize-1;
  warning(['Window size should be even: changing as:',num2str(winsize)]);
 end
 for j=1:numel(num_delay)
  for t=1:(num_delay-winsize-1)
   idx=t:(t+winsize-1);
   ref_time=median(idx*TR_sec);
   postonset_sec=[postonset_sec ref_time];
   cntName=[cntName, [JOB.cntname_fir{k},'[',pad(ref_time,3),'s]']];
   
   cntVec=[];
   for c=1:JOB.num_cond
    x=zeros(1,num_delay);
    x(idx) = ones(1,numel(idx))*JOB.cntvec_fir(k,c);
    cntVec=[cntVec x];
   end
   cntMtx=[cntMtx; cntVec];
  end
 end
end
hf=figure; imagesc(cntMtx);
drawnow;
screen2png([JOB.dir_glm,'/fir_',num2str(length_sec),'s_', ...
 num2str(num_delay),'-th.png']); close(hf);

JOB2=JOB;
JOB2.conname = cntName;
JOB2.convec  = cntMtx;

JOB = myspm_mycon (JOB2);
JOB.postonset_sec = postonset_sec;
end
