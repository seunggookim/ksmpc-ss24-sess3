function [conn,ftchannames] = nbr2connmat(nbr, eeg)

if isfield(eeg,'chanlocs')
  chanlocs = eeg.chanlocs;
else
  chanlocs = eeg;
end

eegchannames = upper({chanlocs.labels}); % including chans wo pos
conn = nan(numel(eegchannames));
ftchannames = upper({nbr.label}); % has channels only with positions
nchans = numel(ftchannames);
assert( sum(ismember(eegchannames,ftchannames)) == numel(eegchannames) )

for i = 1:nchans
  neighbors = upper(nbr(i).neighblabel);
  [~,isnbr] = ismember(neighbors, eegchannames);
  [~,isnnbr] = setdiff(eegchannames, neighbors);
  conn(i,isnbr) = true;
  conn(i,isnnbr) = false;
end
conn = 0.5*(conn+conn');
end