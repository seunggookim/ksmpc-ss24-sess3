function [HO,nii]=myfsl_subcort
nii = load_uns_nii([getenv('FSLDIR'),'/data/atlases/HarvardOxford', ...
 '/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz']);
HO.roi=zeros(size(nii.img));
%HO.roi(ismember(HO.roi, 1+[0 1 2 11 12 13]))=0;
label= 1+[3:6 8:10 14:20 7];
HO.longnames={};
str={'thalamus','caudate','putamen','pallidum','hippocampus','amygdala', 'nucleus accumbens'};
side={'Left','Right'};
k=1;
for s=1:2
for j=1:numel(str)
 HO.longnames{k}=[side{s},' ',str{j}];
 HO.roi(nii.img==label(k)) = k;
 k=k+1;
end
end
HO.longnames=[HO.longnames 'Brainstem'];

str={'Thal','Caud','Ptm','GP','HC','Amyg','NAc'};
side={'-L','-R'};
k=1;
for s=1:2
for j=1:numel(str)
 HO.shortnames{k}=[str{j},side{s}];
 k=k+1;
end
end
HO.shortnames=[HO.shortnames 'Bstem'];

%HOname=xmlread([getenv('FSLDIR'),'/data/atlases/HarvardOxford-Subcortical.xml']);
%  <label index="0" x="58" y="37" z="50">Left Cerebral White Matter</label>
% <label index="1" x="70" y="63" z="35">Left Cerebral Cortex </label>
% <label index="2" x="57" y="40" z="41">Left Lateral Ventrical</label>
% <label index="3" x="51" y="51" z="39">Left Thalamus</label>
% <label index="4" x="51" y="71" z="38">Left Caudate</label>
% <label index="5" x="56" y="67" z="34">Left Putamen</label>
% <label index="6" x="54" y="62" z="35">Left Pallidum</label>
% <label index="7" x="44" y="49" z="18">Brain-Stem</label>
% <label index="8" x="59" y="54" z="27">Left Hippocampus</label>
% <label index="9" x="57" y="61" z="27">Left Amygdala</label>
% <label index="10" x="50" y="70" z="33">Left Accumbens</label>
% <label index="11" x="29" y="38" z="51">Right Cerebral White Matter</label>
% <label index="12" x="25" y="42" z="61">Right Cerebral Cortex </label>
% <label index="13" x="35" y="45" z="44">Right Lateral Ventricle</label>
% <label index="14" x="38" y="51" z="39">Right Thalamus</label>
% <label index="15" x="39" y="72" z="37">Right Caudate</label>
% <label index="16" x="34" y="68" z="34">Right Putamen</label>
% <label index="17" x="35" y="61" z="35">Right Pallidum</label>
% <label index="18" x="31" y="57" z="25">Right Hippocampus</label>
% <label index="19" x="32" y="63" z="25">Right Amygdala</label>
% <label index="20" x="40" y="69" z="32">Right Accumbens</label>

end