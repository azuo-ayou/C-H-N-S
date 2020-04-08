function [Pb,Tb] = genertate_PbTb(P,T,basis_type_trial,N1,N2,x_left,x_right,y_left,y_right)
if basis_type_trial==201
    Pb=P;
    Tb=T;
else if basis_type_trial==202
        % N=2;
        % N1=N;
        % N2=N;
        Pb=zeros(2,(2*N1+1)*(2*N2+1));
        % x_right=1;
        % x_left=0;
        % y_right=1;
        % y_left=0;
        h_x=(x_right-x_left)/(2*N1);
        h_y=(y_right-y_left)/(2*N2);
        
        for k = 1 : 2*N1+1
            Pb(1,(k-1)*(2*N2+1)+1:k*(2*N2+1))=(k-1)*h_x+x_left;
        end
        xx=zeros(2*N1+1,1);
        for i = 1 : 2*N2+1
            for j = 1 : 2*N1+1
                xx(j)=(j-1)*(2*N2+1)+1+(i-1);
            end
            Pb(2,xx)=(i-1)*h_x+y_left;
        end
        Tb=zeros(6,(2*N1*N2));
        B1=zeros(1,2*N2);
        B1(1)=1;
        B1(2*N2)=2*N2+1;
        for i = 2 : 2 : 2*N2-1
            B1(i)=i+1;
            B1(i+1)=i+1;
        end
        for i=1:N1
            B=B1+(i-1)*(2*N2+1)*2;
            for j=1:2*N2
                if mod(j,2)==1;
                    Tb(1,(i-1)*2*N2+j)=B(j);
                    Tb(2,(i-1)*2*N2+j)=B(j)+2*(2*N2+1);
                    Tb(3,(i-1)*2*N2+j)=B(j)+2;
                    Tb(4,(i-1)*2*N2+j)=B(j)+2*N2+1;
                    Tb(5,(i-1)*2*N2+j)=B(j)+2*N2+2;
                    Tb(6,(i-1)*2*N2+j)=B(j)+1;
                else  mod(j,2)==0;
                    Tb(1,(i-1)*2*N2+j)=B(j);
                    Tb(2,(i-1)*2*N2+j)=B(j)+2*(2*N2+1)-2;
                    Tb(3,(i-1)*2*N2+j)=B(j)+2*(2*N2+1);
                    Tb(4,(i-1)*2*N2+j)=B(j)+2*N2;
                    Tb(5,(i-1)*2*N2+j)=B(j)+2*(2*N2+1)-1;
                    Tb(6,(i-1)*2*N2+j)=B(j)+2*N2+1;
                end
            end
        end
    end
end