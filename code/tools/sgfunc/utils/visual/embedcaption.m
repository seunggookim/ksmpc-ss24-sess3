function FnameOut = embedcaption(FnameSrc, Text, Cfg)
if ~exist('Cfg','var'), Cfg= []; end

I = imread(FnameSrc);
Scale = Cfg.WidthPx / size(I,2);
I = imresize(I,Scale);
FnImg = [tempname,'.png'];
imwrite(I, FnImg);
CleanerImg = onCleanup(@() delete(FnImg));
bgcolor = double(squeeze(I(1,1,:)));

hfig = figure('position',[1 1 Cfg.WidthPx 500],'visible','off');
c = uicontrol('Style','Text');
if ~iscell(Text)
  c.String = {Text};
else
  c.String = Text;
end
c.HandleVisibility='off';

if ~isfield(Cfg,'FontSize'), Cfg.FontSize = 12; end
if ~isfield(Cfg,'MarginWidth'), Cfg.MarginWidth = 1*Cfg.FontSize; end
if ~isfield(Cfg,'MarginHeight'), Cfg.MarginHeight = 0.7*Cfg.FontSize; end
c.Position = [Cfg.MarginWidth Cfg.MarginHeight ...
  Cfg.WidthPx-2*Cfg.MarginWidth 500-2*Cfg.MarginHeight];

c.BackgroundColor = bgcolor;
set(gcf,'color',bgcolor);
c.ForegroundColor = 1-bgcolor;
c.FontSize = Cfg.FontSize;
c.HorizontalAlignment = 'Left';

% Fit the UI position to the wrapped text:
[~,pos] = textwrap(c, c.String);
c.Position = pos;
hfig.Position(4) = c.Position(4) + 2*Cfg.MarginHeight;
drawnow

% Get the current DPI:
set(0,'Units','Pixels')
ScreenSize_pixel = get(0,'ScreenSize');
set(0,'Units','Inches')
ScreenSize_inch = get(0,'ScreenSize');
DPI = ScreenSize_pixel ./ ScreenSize_inch;

% save
warning off
FnCap = [tempname,'.png'];
ClenarCap = onCleanup(@() delete(FnCap));
export_fig(FnCap,'-nocrop',sprintf('-r%i',round(DPI(3))))
close(hfig)

%% Now concatenate two images:
[Path,Name,Ext] = fileparts(FnameSrc);
FnameOut = fullfile(Path,[Name,'_cap',Ext]);
imageconcat({FnImg,FnCap}, FnameOut, 1);

warning on
end
