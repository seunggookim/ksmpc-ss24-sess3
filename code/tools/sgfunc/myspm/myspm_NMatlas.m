function [strc,S] = myspm_NMatlas(mni_xyz, isijk)
% [struct_name]=myspm_NMatlas(mni_xyz, isijk)
% mni_xyz: MNI-coordinates (mm) or 1-based ijk (voxels) with nifti
% [1x3] find the name for a given point
% [Nx3] find the most probable name from mode(labels)
%
% with no input arguments, will try to read a coordinate from SPM figure.
% can be used with myspm_graph.m
%
% This function reads labels and returns structure name with the maximal probability
% from [spm12]/tpm/labels_Neuromorphometrics.xml/.nii
% (*) Maximum probability tissue labels derived from the "MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling" (https://masi.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details). These data were released under the Creative Commons Attribution-NonCommercial (CC BY-NC) with no end date. Users should credit the MRI scans as originating from the OASIS project (http://www.oasis-brains.org/) and the labeled data as "provided by Neuromorphometrics, Inc. (http://Neuromorphometrics.com/) under academic subscription".  These references should be included in all workshop and final publications.
%
% (cc) 2015. sgKIM. solleo@gmail.com

%%
if ~exist('isijk','var')
isijk=0;
end

%% read .XML files (from spm12 directory? or somewhere else?)
[spmpath,~,~] = fileparts(which('spm'));
srcpath=[spmpath,'/tpm/'];
xmlFNAME=fullfile(srcpath,'/labels_Neuromorphometrics.xml');
% The problem of the original .XML file is that it contains a comment line 17
% and this needs to be removed to read it using xml2struct.m
% so we just use XML file without those comment lines:
system(['sed -n ''/<!--/!p'' ',xmlFNAME,' > /tmp/labels_Neuromorphometrics.xml']);
s = xml2struct('/tmp/labels_Neuromorphometrics.xml');
S=[];
NL=numel(s.atlas.data.label);
for i=1:NL
S.label(i) = str2double(s.atlas.data.label{i}.index.Text);
S.name{i}  = s.atlas.data.label{i}.name.Text;
end

niiFNAME=fullfile(srcpath,'/labels_Neuromorphometrics.nii');

%% convert xyz coordinate to ijk coordinate
if ~isijk
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
ijk = round(xyz2ijk(xyz, niiFNAME))';
else
ijk = mni_xyz';
xyz = round(ijk2xyz(ijk', niiFNAME))';
end


%% now read probs for XYZ using spm_get_data (very efficient when reading only one voxel)

P = spm_vol(niiFNAME);
label = spm_get_data(P, ijk);

strc.name='null';
strc.prob=100;
if numel(label)==1
idx=find(S.label==label);
if ~isempty(idx)
strc.name = S.name{idx};
strc.prob=100;
end
elseif numel(label)>1
% discard zeros from the labels
label(label==0)=[];
idx = mode(label);
if ~isempty(idx) && ~isnan(idx)
strc.name = S.name{S.label == idx};
strc.prob = round(mean(label==idx)*100);
end
end

strc.mni_xyz = xyz';
strc.ijk1 = ijk';

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
% (cc) sgKIM, 2014. solleo@gmail.com


if ischar(nii)
nii = load_uns_nii(nii);
end

T=[nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];

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

%  ijk = xyz2ijk(xyz, nii)
% converts world-to-voxel coordinates.
%
% Inputs:
%   xyz   [Nx3 vector] is world-coordinates (e.g. MNI-coord)
%   nii   the nii structure read using load_untouch_nii.m (or the filename)
% Output:
%   ijk   [Nx3 vector] is 1-based voxel coordinates for MATLAB
%
% see ijk2xyz.m
% (cc) sgKIM, 2014. solleo@gmail.com

function ijk = xyz2ijk(xyz, nii)

if ischar(nii)
hdr= load_untouch_header_only(nii);
else
hdr = nii.hdr;
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

T=[hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; 0 0 0 1];
ijk = inv(T) * [xyz ones(size(xyz,1),1)]';
ijk(4,:)=[];
ijk = ijk'+1;


end
