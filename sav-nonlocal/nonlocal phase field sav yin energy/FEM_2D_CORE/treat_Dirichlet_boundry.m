function [A,b] =treat_Dirichlet_boundry(A,b,boundrynodes,Dirichlet_fun1,Dirichlet_fun2,Pb,t)
%-1:Dirichlet
number_of_boundrynodes=size(boundrynodes,2);
for k = 1 :number_of_boundrynodes
    if boundrynodes(1,k) == -1   
        i = boundrynodes(2,k);
        A(i,:)=0;
        A(i,i)=1;
%         b(i)=feval(Dirichlet_fun1,Pb(1,i),Pb(2,i),t);%×ø±ê
       b(i) = -1  ;
        A(i+size(Pb,2),:)=0;
        A(i+size(Pb,2),i+size(Pb,2))=1;
%         b(i+size(Pb,2))=feval(Dirichlet_fun2,Pb(1,i),Pb(2,i),t);%×ø±ê    
        b(i+size(Pb,2)) = -1;
    end
end
% A(1+size(Pb,2)*2,:)=0;
% A(1+size(Pb,2)*2,1+size(Pb,2)*2)=1;
% b(1+size(Pb,2)*2)=feval(Dirichlet_fun3,Pb(1,1),Pb(2,1),t);%×ø