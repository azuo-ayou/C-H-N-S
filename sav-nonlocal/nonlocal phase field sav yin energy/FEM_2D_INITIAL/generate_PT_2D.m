function [P,T] = generate_PT_2D(N1,N2,x_left,x_right,y_left,y_right)
h_x=(x_right-x_left)/N1;
h_y=(y_right-y_left)/N2;
P=zeros(2,(N1+1)*(N2+1));
T=zeros(3,2*N1*N2);
for k = 1 : N1+1
    P(1,(k-1)*(N2+1)+1:k*(N2+1))=(k-1)*h_x+x_left;
end
xx=zeros(N1+1,1);
for i = 1 : N2+1
for j = 1 : N1+1
    xx(j)=(j-1)*(N2+1)+1+(i-1);
end
P(2,xx)=(i-1)*h_x+y_left;
end
for k = 1:2*N1*N2
    if k/(2*N2)~=floor(k/(2*N2))
        k1=floor(k/(2*N2))+1;
        k2=k-floor(k/(2*N2))*N2*2;
    else 
        k1=floor(k/(2*N2));
        k2=2*N2;
    end
    if k/2 ~= floor(k/2)
        T(1,k)=find_global_number(k1,floor(k2/2)+1,N2);
        T(2,k)=find_global_number(k1+1,floor(k2/2)+1,N2);
        T(3,k)=find_global_number(k1,floor(k2/2)+2,N2);
    else
        T(1,k)=find_global_number(k1,floor(k2/2)+1,N2);
        T(2,k)=find_global_number(k1+1,floor(k2/2),N2);
        T(3,k)=find_global_number(k1+1,floor(k2/2)+1,N2);
    end
end
        
        