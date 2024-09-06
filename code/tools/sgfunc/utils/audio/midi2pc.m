function Pc = midi2pc(Midi)
%Converts a MIDI structure (from READMIDI) to a pitch-class chord (PC) TABLE structure
% REQUIREMENT: https://github.com/kts/matlab-midi
%
% (cc) 2023-07-03. dr.seunggoo.kim@gmail.com

assert(exist('midiInfo','file'), 'INSTALL & PATH: https://github.com/kts/matlab-midi')
Notes = midiInfo(Midi,0);
Onsets = Notes(:,5);
Offsets = Notes(:,6);
NoteNumbers = Notes(:,3);
%%
TimeOnset_s=[]; TimeOffset_s=[]; PitchSet = {}; PitchClassSet={};
QuantizedOnsets = roudn(Onset/0.89)*0.89;
UniqueQuantizedOnsets = unique(QuantizedOnsets);
for iChord = 1:numel(UniqueQuantizedOnsets)
  Ind = UniqueQuantizedOnsets(iChord) == QuantizedOnsets;
  assert(sum(Ind)==4, '# notes ~= 4')
  PitchSet = [PitchSet; unique(NoteNumbers(Ind))'];
  PitchClassSet = [PitchClassSet; unique(mod(PitchSet{end},12))];
  TimeOnset_s = [TimeOnset_s; mean(Onsets(Ind))];
  TimeOffset_s = [TimeOffset_s; mean(Offsets(Ind))];
end
%%