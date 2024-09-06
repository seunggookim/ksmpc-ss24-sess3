function myfsl_colorbar(LABEL, fname_out)
%%
cmapname={'RY','BC'};
for c=1:2
hf=figure('color','k','position',[2058         184         148          61]);
axes('position',[0.1 0.14 .8 .95/2]);
caxis([0 1])
colormap(sgcolormap(cmapname{c}));
hb1=colorbar('location','north');
set(hb1,'color','w', 'xtick',[0 1], 'xticklabel',{'p=0.05','p=0'});
axis off
title(hb1,LABEL{c},'color','w','fontsize',11);
fname=['/tmp/cb',num2str(c),'.png'];
screen2png(fname,600);
im{c}=imread(fname);
close(hf)
end
I=zeros(366*2,888,3,'uint8');
I(1:366, :, :) = im{1};
I(367:end, :, :) = im{2};
imwrite(I, fname_out, 'png')
end
