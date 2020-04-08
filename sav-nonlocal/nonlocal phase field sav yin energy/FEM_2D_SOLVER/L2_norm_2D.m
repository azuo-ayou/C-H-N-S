function error = L2_norm_2D(solution,exact_solution_fun,Tb,T,P,basis_type,t)
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
    local_error=Guass_int_for_error_2D(Gauss_weight,Gauss_point,exact_solution_fun,uh_local,vertices,basis_type,Tb,t);
    error=error+local_error;
  
end
error=sqrt(error);
end

