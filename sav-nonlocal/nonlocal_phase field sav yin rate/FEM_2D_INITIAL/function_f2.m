function r=function_f2(x,y,t)
ef = 1;

%   r = (256.*x.^2.*y.^2.*cos(pi.*t).*(x - 1).^2.*(y - 1).^2 - 16777216.*x.^6.*y.^6.*cos(pi.*t).^3.*(x - 1).^6.*(y - 1).^6)/ef.^2 + 512.*x.^2.*cos(pi.*t).*(x - 1).^2.*(y - 1).^2 + 512.*y.^2.*cos(pi.*t).*(x - 1).^2.*(y - 1).^2 + 512.*x.^2.*y.^2.*cos(pi.*t).*(x - 1).^2 + 512.*x.^2.*y.^2.*cos(pi.*t).*(y - 1).^2 + 256.*x.^2.*y.^2.*cos(pi.*t).*(x - 1).^2.*(y - 1).^2 + 1024.*x.*y.^2.*cos(pi.*t).*(2.*x - 2).*(y - 1).^2 + 1024.*x.^2.*y.*cos(pi.*t).*(2.*y - 2).*(x - 1).^2;
 
   r = 0;




end