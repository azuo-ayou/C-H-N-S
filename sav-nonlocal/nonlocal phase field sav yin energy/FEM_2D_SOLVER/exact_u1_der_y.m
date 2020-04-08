% function result=exact_u_der_y(x,y)
%  result=(x-0.5*x.*x).*(1-y).*exp(x+y)+((x.*y-0.5*x.*x.*y).*(-y.*exp(x+y)));
% end
function result=exact_u1_der_y(x,y,t)
%    result=-2*pi*sin(2*pi*y).*cos(2*pi*x)*cos(2*pi*t);
%   result=2*x.^2.*y.*cos(2*pi*t)*(x - 1).^2.*(2*y.^2 - 3*y + 1);
 
% result= 512*x.^2.*y.*cos(2*pi*t)*(x - 1).^2.*(2*y.^2 - 3*y + 1);
result= 512.*x.^2.*y.*cos(pi.*t).*(x - 1).^2.*(2.*y.^2 - 3.*y + 1);
  
 %µÚÈý¸ö
% result=-cos(x).*sin(y)*cos(2*pi*t);
end