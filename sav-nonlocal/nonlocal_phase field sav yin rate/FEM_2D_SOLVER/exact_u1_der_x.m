% function result=exact_u_der_x(x,y)
%  result=(y-x.*y).*(1-y).*exp(x+y)+(x.*y-0.5*x.*x.*y).*(1-y).*(exp(x+y));
% end
function result=exact_u1_der_x(x,y,t)
%    result=-2*pi*cos(2*pi*y).*sin(2*pi*x)*cos(2*pi*t);
%  result=2*x.*y.^2.*cos(2*pi*t)*(y - 1).^2.*(2*x.^2 - 3*x + 1);

 
%   result=512*x.*y.^2.*cos(2*pi*t)*(y - 1).^2.*(2.*x.^2 - 3*x + 1);
  result=512.*x.*y.^2.*cos(pi.*t).*(y - 1).^2.*(2.*x.^2 - 3.*x + 1);
 
%µÚÈý¸ö
%  result=-sin(x).*cos(y)*cos(2*pi*t);
end