function error = Hx_norm_2D(solution,exact_der_solution_fun,exact_solution_fun,Tb,T,P,basis_type,t)

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
    local_Hx_error=Guass_int_for_error_Hx_2D(Tb,Gauss_weight,Gauss_point,exact_der_solution_fun,exact_solution_fun,uh_local,vertices,basis_type,t);
    error=error+local_Hx_error;
end
error=sqrt(error);
end

function local_Hx_error=Guass_int_for_error_Hx_2D(Tb,Gauss_weight,Gauss_point,exact_der_solution_fun,exact_solution_fun,uh_local,vertices,basis_type,t)
%����ֲ������ϵ���ֵ����ֵ����Ĳ�ֵ��ƽ���Ļ���ֵ���������Ը�����
%   local_error=Guass_int_for_error_1D(Gauss_weight,Gauss_point,exact_solution_fun,uh_local,vertices,basis_type,der)
%����
int_value=0;
int_value_Hx=0;
Gpn=size(Gauss_point,2);
number_of_local_basis = size(Tb,1);
for i=1:Gpn
    sum=0;
    sumHx=0;
    for j=1:number_of_local_basis
   %int_value=int_value+Gauss_weight(i)*(feval(exact_solution_fun,Gauss_point(i))-FE_function_1D(Gauss_point(i),uh_local,vertices,basis_type,der))^2; 
    sumHx=sumHx+uh_local(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,1,0);
    sum=sum+uh_local(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,0,0);
    end
    int_value_Hx=int_value_Hx+Gauss_weight(i)*(feval(exact_der_solution_fun,Gauss_point(1,i),Gauss_point(2,i),t)-sumHx)^2;
     int_value=int_value+Gauss_weight(i)*(feval(exact_solution_fun,Gauss_point(1,i),Gauss_point(2,i),t)-sum)^2;
    
end
local_Hx_error=int_value+int_value_Hx;
end