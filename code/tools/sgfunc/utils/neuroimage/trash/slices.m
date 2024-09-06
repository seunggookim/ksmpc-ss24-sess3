function h = slices(fn_im, cfg, fn_im2)
% h = slices(fn_image, [cfg], fn_image2)
%
% creates slice montage of a given filename
%
% by default, 3x3 sagittal, coronal, axial slices
%
% cfg may have:
% (.layout)     [1x2]  # of rows (R) and # of columns (C)
% (.slicesdim)  [RxC]  slice dimension [1|2|3]
% (.vcoords)    [RxC]  voxel-coordinates to select slices
% (.caxis)      [1x2]  color range
% (.caxisprc)   [1x2]
% (.z)          [1x1]
% if ~isfield(cfg,'fname_png'), FV='on'; else, FV='off'; end
% if ~isfield(cfg,'dpi'), cfg.dpi=120; end
% if ~isfield(cfg,'layout'), cfg.layout=[3 3]; end
% if ~isfield(cfg,'slicesdim'), cfg.slicesdim=[1 1 1; 2 2 2; 3 3 3]; end
% if ~isfield(cfg,'scale'), cfg.scale=1; end
% if ~isfield(cfg,'textcolor'), cfg.textcolor='w'; end
% if ~isfield(cfg,'figurecolor'), cfg.figurecolor='k'; end
% if ~isfield(cfg,'colormap'), cfg.colormap=gray; end
% if ~isfield(cfg,'caxisprc'), cfg.caxisprc=[1 99]; end
% if ~isfield(cfg,'colorbar'), cfg.colorbar=true; end
% if ~isfield(cfg,'contournum'), cfg.contournum=2; end
% if ~isfield(cfg,'contourcolor'), cfg.contourcolor='r'; end
% if ~isfield(cfg,'linewidth'), cfg.linewidth=0.5; end
% if isfield(cfg,'t'), t=cfg.t; end
%
% (cc) 2018-2019, sgKIM  solleo@gmail.com
%%
if exist('cfg','var') && isnumeric(cfg) && ~isempty(cfg), t=cfg; cfg=[]; end
if ~exist('cfg','var'), cfg=[]; end
if ~isfield(cfg,'fname_png'), FV='on'; else, FV='off'; end
if ~isfield(cfg,'dpi'), cfg.dpi=120; end
if ~isfield(cfg,'layout'), cfg.layout=[3 3]; end
if ~isfield(cfg,'slicesdim'), cfg.slicesdim=[1 1 1; 2 2 2; 3 3 3]; end
if ~isfield(cfg,'scale'), cfg.scale=1; end
if ~isfield(cfg,'textcolor'), cfg.textcolor='w'; end
if ~isfield(cfg,'figurecolor'), cfg.figurecolor='k'; end
if ~isfield(cfg,'colormap'), cfg.colormap = gray; end
if ~isfield(cfg,'caxisprc'), cfg.caxisprc = [1 99]; end
if ~isfield(cfg,'colorbar'), cfg.colorbar = true; end
if ~isfield(cfg,'contournum'), cfg.contournum = 2; end
if ~isfield(cfg,'contourcolor'), cfg.contourcolor = 'r'; end
if ~isfield(cfg,'linewidth'), cfg.linewidth = 0.5; end
if ~isfield(cfg,'tightaxisprc'), cfg.tightaxisprc = 0; end
if isfield(cfg,'t'), t=cfg.t; end
if isfield(cfg,'figurevisible'), FV = cfg.figurevisible; end
if islogical(fn_im)
  im = single(fn_im);
elseif isnumeric(fn_im)
  im = fn_im;
else
  nii = load_untouch_nii(fn_im);
  im = nii.img;
end
d=size(im);
if numel(d)==4
  if ~exist('t','var'), t=1; end
  warning('[!] 4D image given: only showing the %i-th volume',t);
  im=im(:,:,:,t);
end
if ~isfield(cfg,'vcoords')
  cfg.vcoords=zeros(3);
  for j=1:3
    cfg.vcoords(j,:)=round([.4 .5 .6]*d(j));
  end
end
if ~isfield(cfg,'caxis')
  if numel(unique(im(:)))<10
    cfg.caxis = [min(im(:)) max(im(:))];
  else
    cfg.caxis = [prctile(im(:),cfg.caxisprc(1)) ...
      prctile(im(:),cfg.caxisprc(2))];
  end
end
if ~isnan(cfg.layout)
  ij = [median(1:cfg.layout(1)) median(1:cfg.layout(2))];
  colorbarposition = sub2ind(fliplr(cfg.layout), ij(2), ij(1));
  figurepos = [1 1 150*cfg.layout(2)*cfg.scale 150*cfg.layout(1)*cfg.scale];
else
  colorbarposition = nan;
end
%% image2 for contour
if exist('fn_im2','var')
  if isnumeric(fn_im) && isnumeric(fn_im2)
    im2 = fn_im2;
  else
    nii2 = load_untouch_nii_like(fn_im2, fn_im);
    im2 = nii2.img;
  end
end
%%
if ~isfield(cfg,'figurehandle')
  cfg.figurehandle = figure('position',figurepos, 'visible', FV, 'color',cfg.figurecolor);
else
  set(cfg.figurehandle, 'visible', FV, 'color',cfg.figurecolor)
end

if ~isnan(cfg.layout)
  ax = axeslayout(cfg.layout, [0 0 0 0], [0 0 0 0]);
end
S =  reshape(cfg.slicesdim', [], 1);
if numel(cfg.slicesdim) == 1
  cfg.slicesdim = repmat( cfg.slicesdim, [numel(cfg.vcoords) 1] );
end
V = reshape(cfg.vcoords', [], 1);
h = [];
for k = 1:numel(V)
  if ~isnan(cfg.layout)
    h(k) = axespos(ax,k);
  end
  %% bounding box
  if isfield(cfg,'bb')
    if cfg.bb == 1
      IJK = find3(~~im);
      bb = [min(IJK)-1; max(IJK)+1];
    else
      bb = cfg.bb;
    end
    switch cfg.slicesdim(k)
      case 1
        cfg.axis=[bb(:,2)' bb(:,3)'];
      case 2
        cfg.axis=[bb(:,1)' bb(:,3)'];
      case 3
        cfg.axis=[bb(:,1)' bb(:,2)'];
    end
  end
  %%
  switch S(k)
    case 1
      I = squeeze(im(V(k),:,:))';
    case 2
      I = squeeze(im(:,V(k),:))';
    case 3
      I = squeeze(im(:,:,V(k)))';
  end
  imagesc(I);
  if exist('im2','var')
    switch S(k)
      case 1
        I2 = squeeze(im2(V(k),:,:))';
      case 2
        I2 = squeeze(im2(:,V(k),:))';
      case 3
        I2 = squeeze(im2(:,:,V(k)))';
    end
    hold on
    contour(I2, cfg.contournum, 'color',cfg.contourcolor, ...
      'linewidth',cfg.linewidth);
  end
  set(gca,'ydir','nor'); axis off; axis image;
  if isfield(cfg,'axis')
    axis(cfg.axis);
  end
  if cfg.tightaxisprc
    if exist('im2','var')
      axistight(I2, cfg.tightaxisprc);
    else
      axistight(I, cfg.tightaxisprc);
    end
  end
  caxis(cfg.caxis);
  if cfg.colorbar && (k==colorbarposition)
    cb = colorbaro('location','south', 'color',cfg.textcolor);
    pos = cb.Position;
    cb.Position = [pos(1:3) pos(4)*0.4];
    if isfield(cfg,'colorbartitle')
      title(cb,cfg.colorbartitle,'color',cfg.textcolor)
    end
  end
end
colormap(cfg.colormap)
if isfield(cfg,'fname_png')
  export_fig(cfg.fname_png, '-png', '-nocrop', ['-r',num2str(cfg.dpi)]);
  close(cfg.figurehandle)
end
if ~nargout
  clear h
end

end

function axis0 = axistight(slice1, prc)
% axis0 = axistight(slice1)
% (cc) 2019, sgKIM.
slice1 = slice1 > prctile(slice1(:), prc);
axis0 = [
  find(nansum(slice1,1),1,'first')
  find(nansum(slice1,1),1,'last')
  find(nansum(slice1,2),1,'first')
  find(nansum(slice1,2),1,'last')
  ];
axis(axis0);
end

