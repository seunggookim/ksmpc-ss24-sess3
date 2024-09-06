function mp4ani(fname, frames, fps)

if not(exist('fps','var')) || isempty(fps)
  fps = 30;
end

v = VideoWriter(strrep(fname,'.mp4',''), 'MPEG-4');
v.FrameRate = fps;
open(v)
writeVideo(v, frames)
close(v)

end

%%
function test()
figure
Pos = get(gcf, 'Position');
[x,y] = meshgrid(-2:.2:2,-2:.2:2);
z = x .* exp(-x.^2 - y.^2);
[u,v,w] = surfnorm(x,y,z);
quiver3(x,y,z,u,v,w, Color='k'); 
hold on, surf(x,y,z), hold off
set(gca,CameraViewAngleMode='manual', xlim=[-2 2])
frames = [];
for i = linspace(0, 360, 360*10) -120
  view([i, 53])
  drawnow
  frames = [frames getframe(gcf, [0 0 Pos(3) Pos(4)])];
end
end
