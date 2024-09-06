function compare_v1 (fn_pre, fn_post)
% compare_v1 (fn_pre, fn_post)
% (cc) 2019, sgKIM.

[p1,f1,e1] = myfileparts(fn_pre);
fn_pre = [p1,'/',f1,e1];
fa_pre = load_untouch_nii([p1,'/',f1(1:end-2),'FA',e1]);
v1_pre = load_untouch_nii([p1,'/',f1(1:end-2),'V1',e1]);
[p2,f2,e2] = myfileparts(fn_post);
fn_post = [p2,'/',f2,e2];
fa_post = load_untouch_nii([p2,'/',f2(1:end-2),'FA',e2]);
v1_post = load_untouch_nii([p2,'/',f2(1:end-2),'V1',e2]);

fn_gif = [p1,'/',f1,'_vs_',f2,'.gif'];
figure('color','k', 'position',[642 439 1279 253],'visible','off')
Z = round(linspace(1, fa_pre.hdr.dime.dim(4), 19+4));
Z = Z(6:10+2);
Title = {fn_pre, fn_post};
ax = axeslayout([1 7],[0 0 0.05 0],[0 0 0 0]);
I = {v1_pre.img, v1_post.img};
J = {fa_pre.img, fa_post.img};
delay_sec = [2 2];
for i=1:2
  clf
  for j=1:7
    axespos(ax, j)
    rgbslice = permute(abs(squeeze(I{i}(:,:,Z(j),:))),[2 1 3]);
    image(rgbslice .* 1.5 .* zeroone(repmat(J{i}(:,:,Z(j))',[1 1 3])))
    axis off; axis image; set(gca,'ydir','nor')
    axis([
      find(sum(sum(rgbslice,3),1),1,'first')
      find(sum(sum(rgbslice,3),1),1,'last')
      find(sum(sum(rgbslice,3),2),1,'first')
      find(sum(sum(rgbslice,3),2),1,'last')]);
    if j==1
      h = title(Title{i},'color','w', 'interp','none','fontweight','normal', ...
        'horizontalAlignment','left');
      xlim0 = xlim;
      h.Position(1) = xlim0(1);
    end
  end
  gifani(gcf, fn_gif, i, delay_sec(i));
end
end