function result=exact_u2_der_y(x,y,t) 
ef =1;

 
 
 result= - (256.*x.^2.*y.^2.*cos(pi.*t).*(2.*y - 2).*(x - 1).^2 - 100663296.*x.^6.*y.^5.*cos(pi.*t).^3.*(x - 1).^6.*(y - 1).^6 - 100663296.*x.^6.*y.^6.*cos(pi.*t).^3.*(x - 1).^6.*(y - 1).^5 + 512.*x.^2.*y.*cos(pi.*t).*(x - 1).^2.*(y - 1).^2)/ef.^2 - 1536.*x.^2.*cos(pi.*t).*(2.*y - 2).*(x - 1).^2 - 512.*y.^2.*cos(pi.*t).*(2.*y - 2).*(x - 1).^2 - 3072.*x.^2.*y.*cos(pi.*t).*(x - 1).^2 - 1024.*x.^2.*y.*cos(pi.*t).*(y - 1).^2 - 512.*x.^2.*y.^2.*cos(pi.*t).*(2.*y - 2) - 1024.*y.*cos(pi.*t).*(x - 1).^2.*(y - 1).^2 - 2048.*x.*y.*cos(pi.*t).*(2.*x - 2).*(y - 1).^2 - 1024.*x.*y.^2.*cos(pi.*t).*(2.*x - 2).*(2.*y - 2);
 


%result= 512.*x.^2.*y.*cos(pi.*t).*(x - 1).^2.*(2.*y.^2 - 3.*y + 1);


end