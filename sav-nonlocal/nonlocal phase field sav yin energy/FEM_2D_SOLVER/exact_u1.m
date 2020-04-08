function  result=exact_u1(x,y,t,phi)
ef =0.02;
%    result=cos(2*pi*y).*cos(2*pi*x)*cos(2*pi*t);
%    result=x.^2.*(x-1).^2.*y.^2.*(y-1).^2*cos(2*pi*t);
% result = 0*x.*y*t;


%能量

example = 2;

if example==1
    
    %无真解
    result = 0.5*sin(2*pi*x).*sin(2*pi*y)+0.1;
    
    
    
    
    %   r = 4;
    % % 圆变方
    % for i = 1 :size(x,2)
    %
    %         if abs(x(i)-1/2)<1/r && abs(y(i)-1/2)<1/r
    %             result(i) = 1;
    %         elseif abs(x(i)-1/2)<1/r && y(i)==1/2-1/r %abs(y(i)-1/2)==1/r
    %             result(i) = 1;
    %         elseif x(i)==1/2-1/r && abs(y(i)-1/2)<1/r
    %             result(i) = 1;
    %         elseif x(i)==1/2-1/r && y(i)==1/2-1/r
    %             result(i) = 1;
    %
    %         elseif abs(x(i)-1/2)==1/r && abs(y(i)-1/2)<=1/r
    %             result(i) = 0;
    %         elseif abs(x(i)-1/2)<1/r && abs(y(i)-1/2)==1/r
    %             result(i) = 0;
    %         else
    %             result(i) = -1;
    %         end
    % end
    
elseif example==2
    r = 4;
    % 圆变方
    for i = 1 :size(x,2)
        
        if abs(x(i))<1/r && abs(y(i))<1/r
            result(i) = 1;
        elseif abs(x(i))<1/r && y(i)==-1/r %abs(y(i)-1/2)==1/r
            result(i) = 1;
        elseif x(i)==-1/r && abs(y(i))<1/r
            result(i) = 1;
        elseif x(i)==-1/r && y(i)==-1/r
            result(i) = 1;
            
        elseif abs(x(i))==1/r && abs(y(i))<=1/r
            result(i) = 0;
        elseif abs(x(i))<1/r && abs(y(i))==1/r
            result(i) = 0;
        else
            result(i) = -1;
        end
    end
    
    
    
    
    % for i=1:length(x)
    %
    % if 0.3<x(i)&&x(i)<0.7 && 0.3<y(i)&&y(i)<0.7
    %     result(i,1)=1;
    % elseif x(i)==0.3 && 0.3<y(i)&&y(i)<0.7
    %     result(i,1)=0;
    % elseif x(i)==0.7 && 0.3<y(i)&&y(i)<0.7
    %       result(i,1)=0;
    % elseif y(i)==0.3&& 0.3<x(i)&&x(i)<0.7
    %         result(i,1)=0;
    % elseif  y(i)==0.7 &&0.3<x(i)&&x(i)<0.7
    %     result(i,1)=0;
    % else
    %     result(i,1)=-1;
    % end
    %
    %
    % end
    % result = result';
    
    
    
elseif example==3
    % unifrnd (a,b)
    phi = 0.5;
    for i = 1 :size(x,2)
        
        %         result(i)=(2*phi-1)+0.02*unifrnd (-10,10);
        result(i)=0.1*unifrnd (-1,1);
    end
    
    
    
    
    
    
    % result=zeros(1,length(x));
    % for i = 1 :size(x,2)
    %
    %         if 0<x(i)&& x(i)<0.4 && 0.38<y(i)&& y(i)<0.8
    %             a = -tanh((sqrt((x(i)-0.3)^2+(y(i)-0.5)^2)-0.1)/(sqrt(2)*ef));
    %             result(i) = a;
    %         elseif 0.4<x(i)&& x(i)<0.8 && 0.38<y(i)&& y(i)<0.8
    %             a = -tanh((sqrt((x(i)-0.5)^2+(y(i)-0.5)^2)-0.1)/(sqrt(2)*ef));
    %             result(i) = a;
    %
    %         else
    %             a = -tanh((sqrt((x(i)-0.4)^2+(y(i)-0.25)^2)-0.1)/(sqrt(2)*ef));
    %             result(i) = a;
    %         end
    % end
    
    
    
    
elseif example==4
    result=zeros(1,length(x));
    for i = 1 :size(x,2)
        
        if -1<x(i)&& x(i)<0 && -1<y(i)&& y(i)<0.8
            a = -tanh((sqrt((x(i)+0.15)^2+(y(i))^2)-0.15)/(sqrt(2)*ef));
            result(i) = a;
            
        else
            a = -tanh((sqrt((x(i)-0.15)^2+(y(i))^2)-0.15)/(sqrt(2)*ef));
            result(i) = a;
        end
    end
end

end