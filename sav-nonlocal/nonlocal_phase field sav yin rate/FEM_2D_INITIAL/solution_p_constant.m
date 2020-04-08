function error = solution_p_constant(solution,Tb,T,P,basis_type)
error=0;
Gauss_type=9;
number_of_elements=size(T,2);%�������Ԫ����Ŀ
for n=1:number_of_elements
    vertices=P(:,T(:,n));%��ö�Ӧ��Ԫ��ŵĶ�������
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gauss_type);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    
    uh_local=solution(Tb(:,n));%��ö�Ӧ��Ԫ�ľֲ���������ϵ�������õ��ģ�
    %����ֲ������ϵ���ֵ���
    local_error=Guass_int_for_error_2D(Gauss_weight,Gauss_point,uh_local,vertices,basis_type,Tb);
    error=error+local_error;
  
end
error=error/0.25;
end

function local_error=Guass_int_for_error_2D(Gauss_weight,Gauss_point,uh_local,vertices,basis_type,Tb)
%����ֲ������ϵ���ֵ����ֵ����Ĳ�ֵ��ƽ���Ļ���ֵ���������Ը�����
%   local_error=Guass_int_for_error_1D(Gauss_weight,Gauss_point,exact_solution_fun,uh_local,vertices,basis_type,der)
%����$\int_{x_n}^{x_{n+1}} (u^{der}-w_n)^2 dx,w_n=\sum_{k=1}^{N_{lb}} u_{Tb{k,n}}\fai^{der}_{nk}$
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
    int_value=int_value+Gauss_weight(i)*(sum);
end
local_error=int_value;
end
