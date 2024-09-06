function DAT=myft_rmscarpet(DAT)

num_trl=numel(DAT.trial);
d=size(DAT.trial{1});
D=zeros(d(1),d(2),num_trl);
for t=1:num_trl
  D(:,:,t)=DAT.trial{t};
end
DAT.rms=squeeze(rms(D,1));
imagesc(DAT.rms')
colorbar
end