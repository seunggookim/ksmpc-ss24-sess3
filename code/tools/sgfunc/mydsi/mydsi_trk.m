function JOB = mydsi_trk(JOB)
% JOB = mydsi_trk(JOB)
%
%
% ref: http://dsi-studio.labsolver.org/Manual/command-line-for-dsi-studio
%
% (cc) 2017, sgKIM.

outputname=JOB.outputname;

fname_src=JOB.fname_src;
[p1,~,~]=fileparts(fname_src);
if isempty(p1), p1=pwd; end
fname_out=[p1,'/',outputname,'.txt'];
%% Input and output files
cmd=['/Applications/dsi_studio.app/Contents/MacOS/dsi_studio ', ...
 ' --action=trk ', ...
 ' --source=',fname_src, ...
 ' --output=',fname_out];

%% Mask setting:
 %{
SEED: a starting point is randomly selected (at sub-voxel resolution)
      within given seed voxels
ROI:  region of interest. keep tracks that passes an ROI or ALL ROIs
      i.e., inclusive mask, 
ROA:  region of avoidance. discard tracks that passes an ROI or ANY(?) ROIs
      i.e., exclusive mask
END:  keep tracks that "END" (not passing) in this mask
TERM: terminate any tracks that enter this mask
%}
if ~iscell(JOB.roi),  JOB.roi={JOB.roi};   end
num_roi  = numel(JOB.roi);
if ~iscell(JOB.end),  JOB.end={JOB.end};   end
num_end  = numel(JOB.end);
if ~iscell(JOB.term), JOB.term={JOB.term}; end
num_term = numel(JOB.term);

cmd=[cmd,' --seed=',JOB.seed];
cmd=[cmd,' --roi=',JOB.roi{1}];
if num_roi>1
 for i=2:num_roi, cmd=[cmd,' --roi',num2str(i),'=',JOB.roi{i}]; end
end
cmd=[cmd,' --end=',JOB.end{1}]; 
if num_end>1
 for i=2:num_end, cmd=[cmd,' --end',num2str(i),'=',JOB.end{i}]; end
end

%% Tracking control
% defaults:
isRK4=0;          % 0=streamline 1=Fourth-order Runge-Kutta 
seed_count=1e+5;  % # of streamlines
%fa_threshold=0.4;
threshold_index='qa';  % QA (ODF-generalized FA)
qa_threshold=0.04; 
initial_dir=0;    % initial propagation direction 0:primary fiber (default), 1:random, 2:all fiber orientations
min_length=10;    % in mm
max_length=800;
thread_count=4; 

% given options:
if isfield(JOB,'isRK4'), isRK4=JOB.isRK4; end
if isfield(JOB,'seed_count'), seed_count=JOB.seed_count;end
if isfield(JOB,'threshold_index'), threshold_index=JOB.threshold_index; end
if isfield(JOB,'qa_threshold'), qa_threshold=JOB.qa_threshold;end
if isfield(JOB,'initial_dir'), initial_dir=JOB.initial_dir; end
if isfield(JOB,'min_length'), min_length=JOB.min_length; end
if isfield(JOB,'max_length'), max_length=JOB.max_length; end
if isfield(JOB,'thread_count'), thread_count=JOB.thread_count; end


cmd=[cmd, ...
 ' --method=',num2str(isRK4), ...
 ' --seed_count=',num2str(seed_count), ... % or fiber_count ?
 ' --qa_threshold-',num2str(qa_threshold), ...
 ' --initial_dir=',num2str(initial_dir), ...
 ' --min_length=',num2str(min_length), ...
 ' --max_length=',num2str(max_length), ...
 ' --thread_count=',num2str(thread_count), ...
 ];
myunix(cmd);

JOB.fname_out = fname_out;

%{
--seed_plan=0 --interpolation=1 --thread_count=10 
--seed=${dir}/roiText/${seed1} --end=${dir}/roiText/${seed2} 
--seed_count=100000 --nqa_threshold=0.04 --turning_angle=65 
--step_size=.5 --smoothing=0.50 --min_length=10 --max_length=80   
--output=seed1.trk --threshold_index=nqa --fa_threshold=0.04
%}

end