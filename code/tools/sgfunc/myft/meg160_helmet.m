layout=ft_prepare_layout(struct('layout','yokogawa160_helmet.mat'))
figure;
%%
load yokogawa160_helmet.mat lay
clf
scatter(lay.pos(:,1), lay.pos(:,2),'k.')
grid on
axis image
%line(lay.outline{1}(:,1),lay.outline{1}(:,2))
%line(lay.outline{2}(:,1),lay.outline{2}(:,2))
%line(lay.mask{1}(:,1), lay.mask{1}(:,2),'color','b')

theta=[137:1:403]/360*2*pi;
x=cos(theta)*0.88;
y=sin(theta)*0.88+0.03;
line(x,y,'color','b')
x2=[x(end) 0.4 0 -0.4 x(1)];
y2=[y(end) 0.5 0.56 0.5 y(1)];
line(x2,y2,'color','g')

helmet=[[x x2]' [y y2]'];

theta=[127:-1:53]/360*2*pi;
x=cos(theta)*0.63;
y=sin(theta)*0.63+0.00;
% theta=[130:-1:50]/360*2*pi;
% x=cos(theta)*0.68;
% y=sin(theta)*0.68+0.00;
line(x,y)
face=[x' y'];

i=find(x>0,1)
x=[x(i-8) 0 x(i+8)];
y=[y(i-8) 0.70 y(i+8)];
line(x,y)
nose=[x' y'];

lay.outline{1}=face;
lay.outline{2}=nose;
lay.mask{1}=helmet;

save /home/sgk/matlab/fieldtrip-20180903/template/layout/yokogawa160_helmet.mat lay