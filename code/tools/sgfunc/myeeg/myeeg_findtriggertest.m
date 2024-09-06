function [eeg, isTriggerTest] = myeeg_findtriggertest(eeg)
% Finds trigger test (16 triggers from S1 to R1) events and rename them as 
%   'test'.
%
% USAGE
% [eeg, isTriggerTest] = myeeg_findtriggertest(eeg)
%
% (cc) 2020, sgKIM.

testcodes = cell(1,16);
for i = 1:15
  testcodes{i} = sprintf('S% 3i',i);
end
testcodes{16} = 'R  1';
ntrig = numel(testcodes);

eventtypes = {eeg.event.type};
isTriggerTest = false(size(eventtypes));
for itrial = ntrig:numel(isTriggerTest)
  idx = itrial-(ntrig-1):itrial;
  if all(strcmp(eventtypes(idx), testcodes))
    isTriggerTest(idx) = true;
  end
end
if sum(isTriggerTest)
  idx = find(isTriggerTest);
  for i = idx
    eeg.event(i).type = 'test';
  end
end
end