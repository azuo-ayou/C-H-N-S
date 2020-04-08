function local_error=Guass_int_for_error_2D(Gauss_weight,Gauss_point,exact_solution_fun,uh_local,vertices,basis_type,Tb,t)
%计算局部网格上的真值与数值结果的差值的平方的积分值（导数可以给定）
%   local_error=Guass_int_for_error_1D(Gauss_weight,Gauss_point,exact_solution_fun,uh_local,vertices,basis_type,der)
%计算$\int_{x_n}^{x_{n+1}} (u^{der}-w_n)^2 dx,w_n=\sum_{k=1}^{N_{lb}} u_{Tb{k,n}}\fai^{der}_{nk}$
int_value=0;
%Gpn=size(Gauss_point,2);
Gpn=9;
number_of_local_basis = size(Tb,1);
for i=1:Gpn
    sum=0;
    for j=1:number_of_local_basis
   %int_value=int_value+Gauss_weight(i)*(feval(exact_solution_fun,Gauss_point(i))-FE_function_1D(Gauss_point(i),uh_local,vertices,basis_type,der))^2; 
   
    sum=sum+uh_local(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,0,0);
    end
    int_value=int_value+Gauss_weight(i)*(feval(exact_solution_fun,Gauss_point(1,i),Gauss_point(2,i),t)-sum)^2;
end
local_error=int_value;
