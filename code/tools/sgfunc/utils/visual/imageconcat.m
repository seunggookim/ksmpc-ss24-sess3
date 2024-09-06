function Iout= imageconcat(inputs, fname_output, dim, cfg)
% SYNTAX:
% imageconcat(fnames_input, [fname_output], [dim], [cfg])
%
% INPUT:
%  inputs       {1xN} with filenames (char) or images (uint8)
%  fname_output '1xN' output filename
%  dim          [1x1] dimension to concatanate: (1=rows, 2=columns)
%
%  cfg may include:
%   .crop   [1x4] [top bottom left right] in pixels
%   .bgval  [1x1] (255 for white[default], 0 for black)
%   .margin [1x1] space between images in pixels
%
% OUTPUT:
%  Iout is a image matrix
%  if fname_output is given, Iout will be saved
%  if no output argument is defined, it will show Iout for 2 seconds and
%  close it.
%
% (cc) 2017, 2019, sgKIM

%% parse options
if ~exist('dim','var'), dim=1; end
if ~exist('cfg','var'), cfg=[]; end
% if ~isfield(cfg,'isPadding'), cfg.isPadding=1; end
if ~isfield(cfg,'crop'), cfg.crop=[0 0 0 0]; end
if ~isfield(cfg,'bgval'), cfg.bgval=255; end
if ~isfield(cfg,'margin'), cfg.margin=0; end


%% read inputs
N=numel(inputs);
D=zeros(N,2);
Iin=cell(1,N);
for i=1:N
  if isnumeric(inputs{i}) % if numbers
    Iin{i}=inputs{i};
  elseif ischar(inputs{i})
    Iin{i}=imread(inputs{i}); % if filenames
  else
    error('"inputs" cell needs to contain images or filenames');
  end
  % crop
  Iin{i}(1:cfg.crop(1),:,:)=[];
  Iin{i}(end-cfg.crop(2)+1:end ,:,:)=[];
  Iin{i}(:,1:cfg.crop(3),:)=[];
  Iin{i}(:,end-cfg.crop(4)+1:end,:)=[];
  % read dimension
  for d=1:2
    D(i,d)=size(Iin{i},d);
  end
  % if the input image is BW, make it RGB:
  if size(Iin{i},3) == 1
    Iin{i} = repmat(Iin{i},[1 1 3]);
  end
end


%% concatanate images for each channel
dim2match=2*(dim==1) + 1*(dim==2);
dim_max = max(D(:,dim2match));
RGB=cell(1,3);
for i=1:N
  for rgb=1:3
    i1=Iin{i}(:,:,rgb);
    dim1=size(i1,dim2match);
    if dim1 ~= dim_max % need to do something?
      if isfield(cfg,'isPaddingCenter') && cfg.isPaddingCenter % is it zero-padding?
        numpix1=floor((dim_max-dim1)/2);
        numpix2=(dim_max-dim1)-numpix1;
        if dim==1
          z1=cfg.bgval*ones(size(i1,1),numpix1,'uint8');
          z2=cfg.bgval*ones(size(i1,1),numpix2,'uint8');
          i1=[z1, i1, z2];
        else
          z1=cfg.bgval*ones(numpix1,size(i1,2),'uint8');
          z2=cfg.bgval*ones(numpix2,size(i1,2),'uint8');
          i1=[z1; i1; z2];
        end
      elseif isfield(cfg,'isPaddingLeft') && cfg.isPaddingLeft % is it zero-padding?
        numpix1=dim_max-dim1;
        if dim==1
          z1=cfg.bgval*ones(size(i1,1),numpix1,'uint8');
          i1=[i1, z1];
        else
          z1=cfg.bgval*ones(numpix1,size(i1,2),'uint8');
          i1=[i1; z1];
        end
      else % or resizing?
        if dim==1
          i1=imresize(i1,[nan dim_max],'bicubic');
        else
          i1=imresize(i1,[dim_max nan],'bicubic');
        end
      end
    end
    
    if dim == 1
      space = ones(cfg.margin*(i>1), size(RGB{rgb},2))*cfg.bgval;
      RGB{rgb} = [RGB{rgb}; space; i1];
    else
      space = ones(size(RGB{rgb},1), cfg.margin*(i>1))*cfg.bgval;
      RGB{rgb} = [RGB{rgb}, space, i1];
    end
  end
end


%% combine channels
Iout = zeros([size(RGB{1}) 3],'uint8');
for r = 1:3
  Iout(:,:,r) = RGB{r};
end


%% SCALE:
if isfield(cfg,'scale')
  Iout = imresize(Iout, cfg.scale, 'cubic');
end

  
%% OUTPUT:
if exist('fname_output','var')
  imwrite(Iout, fname_output);
end
if nargout == 0
  hf = figure('visible','off'); image(Iout); axis image tight; 
  drawnow; pause(2); close(hf); clear Iout
end
end

