%% Œ»Ã¨
% syms x y niu 
% u1=x^2*y^2+exp(-y);
% u2=(-2/3)*x*y^3+2-pi*sin(pi*x);
% p=(pi*sin(pi*x)-2)*cos(2*pi*y);
% 
% f11=-niu*(diff(u1,x,2)+diff(u1,y,2))+diff(p,x);
% f12=u1*diff(u1,x)+u2*diff(u1,y);
% 
% simplify(f11+f12)
% 
% f21=-niu*(diff(u2,x,2)+diff(u2,y,2))+diff(p,y);
% f22=u1*diff(u2,x)+u2*diff(u2,y);
% simplify(f21+f22)
%% ∑«Œ»Ã¨
% syms x y t niu
% u1=(x^2*y^2+exp(-y))*cos(2*pi*t);
% u2=((-2/3)*x*y^3+2-pi*sin(pi*x))*cos(2*pi*t);
% p=(pi*sin(pi*x)-2)*cos(2*pi*y)*cos(2*pi*t);
% 
% f11=diff(u1,t)-niu*(diff(u1,x,2)+diff(u1,y,2))+diff(p,x);
% f12=u1*diff(u1,x)+u2*diff(u1,y);
% simplify(f11)
% 
% f21=diff(u2,t)-niu*(diff(u2,x,2)+diff(u2,y,2))+diff(p,y);
% f22=u1*diff(u2,x)+u2*diff(u2,y);
% simplify(f21)


%% phase field
syms x y t ef
% u1=cos(2*pi*y).*cos(2*pi*x)*cos(2*pi*t);

 u1 =16*16*x^2*(x-1)^2*y^2*(y-1)^2*cos(pi*t);
%  niu=16*16*x^2*(x-1)^2*y^2*(y-1)^2*cos(pi*t);
% u2=((-2/3)*x*y^3+2-pi*sin(pi*x))*cos(2*pi*t);
% p=(pi*sin(pi*x)-2)*cos(2*pi*y)*cos(2*pi*t);

niu = -(diff(u1,x,2)+diff(u1,y,2))+1/ef^2*(u1^3-u1);
 f1 = diff(u1,t)-(diff(niu,x,2)+diff(niu,y,2));
% f2 = niu+(diff(u1,x,2)+diff(u1,y,2))-1/ef^2*(u1^3-u1);


%  simplify(diff(u1,x))
%  simplify(diff(u1,y))
 simplify(niu)



% syms x y t niu ef
% u1=cos(x)*cos(y)*cos(2*pi*t);
% niu = -(diff(u1,x,2)+diff(u1,y,2))+1/ef^2*(u1^3-u1);
% f12=diff(u1,t)-(diff(niu,x,2)+diff(niu,y,2));
% % simplify(niu)
%  simplify(f12)
% %  simplify(diff(niu,x))
% %  simplify(diff(niu,y))


