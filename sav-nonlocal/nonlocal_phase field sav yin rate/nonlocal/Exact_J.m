function r = Exact_J(x,y)
% alph = 0.1;
% thta = 0.001;%N=64��epsion=0.02�����Բ
% thta = sqrt(0.1);
 thta = sqrt(0.1);

alph = 4/(pi*thta^2);

%   r=alph/thta^2*exp(-((x-0.4).^2+(y-0.4).^2)/thta^2);
 r=alph/thta^2.5*exp(-(x.^2+y.^2)/thta^2);
end

%thta^4���б������,thta^2û�б仯;