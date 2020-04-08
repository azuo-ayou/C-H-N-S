function [A,b] =treat_Dirichlet_boundry(N1,A,b,boundrynodes,Dirichlet_fun,Pb)
%-1:Dirichlet
number_of_boundrynodes=size(boundrynodes,2);
for k = 1 :number_of_boundrynodes
    if boundrynodes(1,k) == -1   
        i = boundrynodes(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i)=feval(Dirichlet_fun,Pb(1,i),Pb(2,i));%зјБъ
        
        A(i+(N1+1)^2,:)=0;
        A(i+(N1+1)^2,i+(N1+1)^2)=1;
        b(i+(N1+1)^2)=feval(Dirichlet_fun,Pb(1,i),Pb(2,i));%зјБъ
        
    end
end