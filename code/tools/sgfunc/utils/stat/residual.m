function [Yres,bhat,Yadj] = residual(Y,M)
% [Yres,bhat] = residual(Y,M)
%
% for vector y (Nx1), M is (NxK) for K parameters
% for matrix y (NxP) for P measures, M is (NxPxK).
%
% It returns:
% Yres = Y-M*((M'*M)\M'*Y);
%
% (cc) 2017, 2019. sgKIM.
% yres = y - M*((M'*M)\M'*y);

% sanity check:
if size(Y,1) ~= size(M,1)
  error('Rows of Y and M should be equal!')
end
if ndims(M) == 3
  if size(Y,2) ~= size(M,2)
    error('For matrix Y, the number of measures should be equal!')
  end
end

meanY = mean(Y,1);
Y = Y - meanY;
M = M - mean(M,1);
if isvector(Y)
  bhat = (M'*M)\M'*Y;
  Yres = Y - M*bhat;
else
  bhat = zeros(size(M,2), size(M,3)); % <P-by-K>
  Yres = Y;
  for p = 1:size(M,2) % for each p measure
    Mi = squeeze(M(:,p,:));
    bhat(p,:) = (Mi'*Mi)\Mi'*Y(:,p);
    Yres(:,p) = Y(:,p) - Mi*bhat(p,:)';
  end
end
Yadj = Yres + meanY;

end


%{
%--------------------------------------------
%TESTIING CODE for 1 measure & 1 variables
clear; clc
x = demean(rand(1000,1));
M = demean(rand(1000,1));
b = 0.5;
y = x + b*M + 3;

clf;
subplot(121)
hold on
scatter(x,x+3,'o')
scatter(x,y)
meany = mean(y);
dey = y - meany;
deM = demean(M);
bhat = inv(deM'*deM)*deM'*dey;
yres = dey - M*bhat;
yres = yres + meany;
% plot(x,betahat*M+100)
scatter(x,yres,'r.')
title(sprintf('truebeta=%.3f betahat=%.3f',b, bhat))

legend({'x','y','yres'},'box','on')

subplot(122)
hold on
scatter(x,x+3,'o')
scatter(x,y)
[yres,bhat] = residual(y, deM);
scatter(x,yres,'r.')
title(sprintf('truebeta=%.3f betahat=%.3f',b, bhat))

legend({'x','y','yres'},'box','on')

%--------------------------------------------
%TESTIING CODE for 3 measure & 1 variables
clear; clc
x = demean(rand(1000,3));
M = demean(rand(1000,3,1));
b = [0.5 0.3 1];
y = [];
for i = 1:3
  y(:,i) = x(:,i) + M(:,i,1)*b(i) + 3;
end
[yres,bhat] = residual(y,M);
figure
for i = 1:3
  subplot(1,3,i); hold on
  scatter(x(:,i),x(:,i)+3,'o')
  scatter(x(:,i),y(:,i))
  scatter(x(:,i),yres(:,i),'r.')  
  title(sprintf('truebeta=%.3f betahat=%.3f',b(i), bhat(i)))
  legend({'x','y','yres'},'box','on')
end

%--------------------------------------------
%TESTIING CODE for 3 measure & 2 variables
clear; clc
x = demean(rand(1000,3));
M = demean(rand(1000,3,2));
b = [0.5 0.3 1; 2 4 6]';
y = [];
for p = 1:3
  y(:,p) = x(:,p) + M(:,p,1)*b(p,1) + M(:,p,2)*b(p,2) + 3;
end
[yres,bhat] = residual(y,M);
figure
for p = 1:3
  subplot(1,3,p); hold on
  scatter(x(:,p),x(:,p)+3,'o')
  scatter(x(:,p),y(:,p))
  scatter(x(:,p),yres(:,p),'r.')  
  title(sprintf('b1=%.3f bhat1=%.3f \n b2=%.3f bhat2=%.3f',...
    b(p,1), bhat(p,1), b(p,2), bhat(p,2) ))
  legend({'x','y','yres'},'box','on')
end


%}