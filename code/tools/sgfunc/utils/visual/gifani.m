function gifani(hf, filename, frameid, delay_sec)
% gifani(figurehandle, filename, frameid, delay_sec)
%
%[Example]
%
% delay_sec=0.1;
% filename='test.gif';
% hf=figure;
% for frameid=1:100
%  gifani(hf, filename, frameid, delay_sec)
% end
% (cc) 2015, sgKIM.   solleo@gmail.com

if ~exist('delay_sec','var')
  delay_sec=1;
end

% figure(hf);
drawnow;
frame = getframe(hf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if frameid==1
  imwrite(imind, cm, filename, 'gif', 'loopcount', inf, 'DelayTime', delay_sec);
else
  imwrite(imind, cm, filename, 'gif', 'writemode', 'append', 'DelayTime', delay_sec);
end
end

%{
[Example]

delay_sec=0.1;
filename='test.gif';
hf=figure;
for frameid=1:100
 gifani(hf, filename, frameid, delay_sec)
end

%}