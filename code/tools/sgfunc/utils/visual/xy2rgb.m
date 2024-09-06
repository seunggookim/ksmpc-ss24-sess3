function rgb = xy2rgb(x,y)
% rgb = xy2rgb(x,y)
%
% converting a 2-D coordinate (x,y) into a RGB value.
%
% the length of the vector (sqrt(x^2+y^2)) should be <=1
%
% (cc) 2022, dr.seunggoo.kim@gmail.com

theta = atan2(y,x);
rho = sqrt(x.^2 + y.^2);
rgb = hsv2rgb([mod(theta,2*pi)/(2*pi), rho, ones(size(rho,1),1)]);


end