function [strc, strc_all]=myfsl_atlasquery(mni_xyz, ijkflag, atlasset)
% [strc, strc_all]=myfsl_atlasquery(mni_xyz, ijkflag)
% mni_xyz: MNI-coordinates (mm) or 1-based ijk (voxels) with nifti
%
% will find name for a given MNI coordiante (x,y,z) mm
%   or voxel index (i,j,k) by marking ijkflag
% with no input arguments, will try to read a coordinate from SPM figure.
%
% (cc) 2015. sgKIM. solleo@gmail.com
%
% SEE ALSO: MYFSL_FINDLABEL

if ~exist('ijkflag','var')
  ijkflag=0;
end
if ~exist('atlasset','var')
  atlasset='gm';
end

fslpath = getenv('FSLDIR');
% xmlFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford-Cortical-Lateralized.xml');
% niiFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford/HarvardOxford-cortl-prob-1mm.nii.gz');
xmlFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford-Cortical.xml');
niiFNAMES{1}=fullfile(fslpath,'data/atlases/HarvardOxford/HarvardOxford-cort-prob-1mm.nii.gz');
xmlFNAMES{2}=fullfile(fslpath,'data/atlases/HarvardOxford-Subcortical.xml');
niiFNAMES{2}=fullfile(fslpath,'data/atlases/HarvardOxford/HarvardOxford-sub-prob-1mm.nii.gz');
xmlFNAMES{3}=fullfile(fslpath,'data/atlases/Cerebellum_MNIfnirt.xml');
niiFNAMES{3}=fullfile(fslpath,'data/atlases/Cerebellum/Cerebellum-MNIfnirt-prob-1mm.nii.gz');
xmlFNAMES{4}=fullfile(fslpath,'data/atlases/JHU-tracts.xml');
niiFNAMES{4}=fullfile(fslpath,'data/atlases/JHU/JHU-ICBM-tracts-prob-1mm.nii.gz');
xmlFNAMES{5}=fullfile(fslpath,'data/atlases/JHU-labels.xml');
niiFNAMES{5}=fullfile(fslpath,'data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz');

% okay unzipping takes more than 2 sec..., so try to find it from /tmp first
for a=1:5
  [~,fname1,ext1]=fileparts(niiFNAMES{a});
  if ~exist(['/tmp/',fname1], 'file')
    gunzip(niiFNAMES{a}, '/tmp/');
  end
  niiFNAMES{a} = ['/tmp/',fname1];
  xDOC{a}=xml2struct(xmlFNAMES{a});
end

%% get ijk coordinate
if ~ijkflag
  if nargin == 0
    hReg= evalin('base','hReg;');
    xSPM= evalin('base','xSPM;');
    if numel(hReg) == 1
      xyz = spm_XYZreg('GetCoords',hReg);
    else
      xyz = hReg;
    end
  else
    xyz = mni_xyz';
  end
  try xSPM.XYZmm
    [xyz,i] = spm_XYZreg('NearestXYZ', xyz ,xSPM.XYZmm);
  catch ME
    xyz = mni_xyz';
  end
  ijk = round(xyz2ijk(xyz, niiFNAMES{1}))';
else
  ijk = mni_xyz';
  xyz = round(ijk2xyz(ijk', niiFNAMES{1}))';
end


%% now read probs for XYZ using spm_get_data (very efficient when reading only one voxel)
strc.name='N/A';
strc.prob=0;
strc_all.name=cell(1,5);
strc_all.prob=[0 0 0 0 0];
k=1;
if strcmpi('all',atlasset)
  ATLAS=1:5;
elseif strcmpi('gm',atlasset)
  ATLAS=1:3;
elseif strcmpi('wm',atlasset)
  ATLAS=4:5;
end

for a=ATLAS % for each atlas
  P = spm_vol(niiFNAMES{a});
  probs = spm_get_data(P, ijk);
  if a==2
    % ignore cerebral grey/white matter in the subcortical atlas
    probs([1 2 12 13],:)=0;
  end
  
  % find maximal prob from
  if size(ijk,2) == 1 % this works for probs<Nx1>
    [~,b]=max(probs);
    probs_b = probs(b);
  else % in case of cluster.. probs<NxV>
    % I want to compute mean prob for each label
    probs_ri = round(mean(probs,2));
    [~,b]=max(probs_ri);
    probs_b = probs_ri(b);
  end
  
  if probs_b > strc.prob % is the current maximal greater than the previous one?
    strc.name = xDOC{a}.atlas.data.label{b}.Text;
    strc.prob = probs_b;
  end
  
  % for any other possibilities...
  
  if size(ijk,2) == 1 % for a voxel
    nz = find(~~probs); % non-zero probs.
    for j=1:numel(nz)
      strc_all.name{k}=xDOC{a}.atlas.data.label{nz(j)}.Text;
      strc_all.prob(k)=probs(nz(j));
      k=k+1;
    end
  else % for a cluster
    nz = find(~~probs_ri);
    for j=1:numel(nz)
      strc_all.name{k}=xDOC{a}.atlas.data.label{nz(j)}.Text;
      strc_all.prob(k)=probs_ri(nz(j));
      k=k+1;
    end
  end
end

strc_unsort=strc_all;
% and sort;
[~,idx] = sort(strc_unsort.prob, 'descend');
for j=1:numel(idx)
  strc_all.name{j} = strc_unsort.name{idx(j)};
  strc_all.prob(j) = strc_unsort.prob(idx(j));
end

end


%% =============================================================================
% source: http://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

if (nargin < 1)
  clc;
  help xml2struct
  return
end

if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
  % input is a java xml object
  xDoc = file;
else
  %check for existance
  if (exist(file,'file') == 0)
    %Perhaps the xml extension was omitted from the file name. Add the
    %extension and try again.
    if (isempty(strfind(file,'.xml')))
      file = [file '.xml'];
    end
    
    if (exist(file,'file') == 0)
      error(['The file ' file ' could not be found']);
    end
  end
  %read the xml file
  xDoc = xmlread(file);
end

%parse xDoc into a MATLAB structure
s = parseChildNodes(xDoc);

end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
% Recurse over node children.
children = struct;
ptext = struct; textflag = 'Text';
if hasChildNodes(theNode)
  childNodes = getChildNodes(theNode);
  numChildNodes = getLength(childNodes);
  
  for count = 1:numChildNodes
    theChild = item(childNodes,count-1);
    [text,name,attr,childs,textflag] = getNodeData(theChild);
    
    if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
      %XML allows the same elements to be defined multiple times,
      %put each in a different cell
      if (isfield(children,name))
        if (~iscell(children.(name)))
          %put existsing element into cell format
          children.(name) = {children.(name)};
        end
        index = length(children.(name))+1;
        %add new element
        children.(name){index} = childs;
        if(~isempty(fieldnames(text)))
          children.(name){index} = text;
        end
        if(~isempty(attr))
          children.(name){index}.('Attributes') = attr;
        end
      else
        %add previously unknown (new) element to the structure
        children.(name) = childs;
        if(~isempty(text) && ~isempty(fieldnames(text)))
          children.(name) = text;
        end
        if(~isempty(attr))
          children.(name).('Attributes') = attr;
        end
      end
    else
      ptextflag = 'Text';
      if (strcmp(name, '#cdata_dash_section'))
        ptextflag = 'CDATA';
      elseif (strcmp(name, '#comment'))
        ptextflag = 'Comment';
      end
      
      %this is the text in an element (i.e., the parentNode)
      if (~isempty(regexprep(text.(textflag),'[\s]*','')))
        if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
          ptext.(ptextflag) = text.(textflag);
        else
          %what to do when element data is as follows:
          %<element>Text <!--Comment--> More text</element>
          
          %put the text in different cells:
          % if (~iscell(ptext)) ptext = {ptext}; end
          % ptext{length(ptext)+1} = text;
          
          %just append the text
          ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
        end
      end
    end
    
  end
end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
% Create structure of node info.

%make sure name is allowed as structure name
name = toCharArray(getNodeName(theNode))';
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');

attr = parseAttributes(theNode);
if (isempty(fieldnames(attr)))
  attr = [];
end

%parse child nodes
[childs,text,textflag] = parseChildNodes(theNode);

if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
  %get the data of any childless nodes
  % faster than if any(strcmp(methods(theNode), 'getData'))
  % no need to try-catch (?)
  % faster than text = char(getData(theNode));
  text.(textflag) = toCharArray(getTextContent(theNode))';
end

end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = struct;
if hasAttributes(theNode)
  theAttributes = getAttributes(theNode);
  numAttributes = getLength(theAttributes);
  
  for count = 1:numAttributes
    %attrib = item(theAttributes,count-1);
    %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
    %attributes.(attr_name) = char(getValue(attrib));
    
    %Suggestion of Adrian Wanner
    str = toCharArray(toString(item(theAttributes,count-1)))';
    k = strfind(str,'=');
    attr_name = str(1:(k(1)-1));
    attr_name = strrep(attr_name, '-', '_dash_');
    attr_name = strrep(attr_name, ':', '_colon_');
    attr_name = strrep(attr_name, '.', '_dot_');
    attributes.(attr_name) = str((k(1)+2):(end-1));
  end
end
end

%% ------ SG's functions

function xyz = ijk2xyz(ijk, nii)
% xyz = ijk2xyz(ijk, nii)
% converts world-to-voxel coordinates.
%
% Inputs:
%   ijk   [Nx3 vector] is 1-based voxel coordinates for MATLAB
% Output:
%   xyz   [Nx3 vector] is world-coordinates (e.g. MNI-coord)
%   nii   the nii structure read using load_untouch_nii.m (or the filename)
%
% see xyz2ijk.m
% (cc) sgKIM, 2014, 2022. drseunggookim@gmail.com


if ischar(nii)
%   nii = load_uns_nii(nii);
  info = niftiinfo(nii);
end

% T=[nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];
T = info.Transform.T';

if numel(ijk) == 3
  if size(ijk,1) > size(ijk,2)
    ijk=ijk';
  end
else
  if size(ijk,1) < size(ijk,2)
    ijk=ijk';
  end
end

xyz = (T)*[(ijk-1) ones(size(ijk,1),1) ]';
xyz(4,:)=[];
xyz=xyz';

end

%%  ijk = xyz2ijk(xyz, nii)
% converts world-to-voxel coordinates.
%
% Inputs:
%   xyz   [Nx3 vector] is world-coordinates (e.g. MNI-coord)
%   nii   the nii structure read using load_untouch_nii.m (or the filename)
% Output:
%   ijk   [Nx3 vector] is 1-based voxel coordinates for MATLAB
%
% see ijk2xyz.m
% (cc) sgKIM, 2014, 2022. drseunggookim@gmail.com

function ijk = xyz2ijk(xyz, nii)

if ischar(nii)
%   hdr= load_untouch_header_only(nii);
  info = niftiinfo(nii);
elseif isstruct(nii)
%   hdr = nii.hdr;
  info = nii.info;
end

if numel(xyz) == 3
  if size(xyz,1) > size(xyz,2)
    xyz=xyz';
  end
else
  if size(xyz,1) < size(xyz,2)
    xyz=xyz';
  end
end

% T=[hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; 0 0 0 1];
T = info.Transform.T';
ijk = inv(T) * [xyz ones(size(xyz,1),1)]';
T = info.Transform.T';

ijk(4,:)=[];
ijk = ijk'+1;

end
