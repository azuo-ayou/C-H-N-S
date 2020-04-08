clear all
% xx=-101:2:101;
% yy=-101:2:101;
% [x,y]=meshgrid(xx,yy);
% figure;
% r1=p_1(x,y);
% % r1(x^2+y^2-1<1)=nan;
% surf(x,y,r1)


% 
% [x,y] = meshgrid(-10:.01:10);
% z = ones(size(x));
% z(sin(x)+y.*sin(y)<1) = nan;
% figure
% mesh(x,y,z)
% view(2)

[x,y,z]=meshgrid(linspace(-10,10,20));
f = x.^2+y.^2-z.^2;
figure
% isosurface(x,y,z,f,0,x);
 fv=isosurface(x,y,z,f,0,x);
 patch(fv);
% hold on
% isosurface(x,y,z,f,1,y);
% hold on
% isosurface(x,y,z,f,1,z);