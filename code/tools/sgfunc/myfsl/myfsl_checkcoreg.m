function JOB = myfsl_checkcoreg (JOB)
% JOB = myfsl_checkcoreg (JOB)
%
%
% (cc) 2017, sgKIM, solleo@gmail.com

% I am using FSL and Freesurfer for this....
fsldir=getenv('FSLDIR');
if isempty(fsldir)
error('This script requires FSL, a correct path needs to be set in $FSLDIR');
end


end
