function error = L2_norm_2D_interpolation_error_int(solution,exact_solution,Tb,T,P,basis_type,t,Nx_max,Ny_max,Tb1,T1,P1,x_left,y_left)
error=0;
Gauss_type=9;
number_of_elements=size(T,2);%�������Ԫ����Ŀ
for n=1:number_of_elements
    vertices=P(:,T(:,n));%��ö�Ӧ��Ԫ��ŵĶ�������
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gauss_type);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    n;
    Gauss_point=Gauss_point';
    
    uh_local=solution(Tb(:,n));%��ö�Ӧ��Ԫ�ľֲ���������ϵ�������õ��ģ�
    %����ֲ������ϵ���ֵ���
    local_error=Guass_int_for_error_2D_interpolation_error_int(Gauss_weight,Gauss_point,exact_solution,uh_local,vertices,basis_type,Tb,t,Nx_max,Ny_max,Tb1,T1,P1,x_left,y_left);
    error=error+local_error;
    
end
error=sqrt(error);
end


function local_error=Guass_int_for_error_2D_interpolation_error_int(Gauss_weight,Gauss_point,exact_solution,uh_local,vertices,basis_type,Tb,t,Nx_max,Ny_max,Tb1,T1,P1,x_left,y_left)
%����ֲ������ϵ���ֵ����ֵ����Ĳ�ֵ��ƽ���Ļ���ֵ���������Ը�����
%   local_error=Guass_int_for_error_1D(Gauss_weight,Gauss_point,exact_solution_fun,uh_local,vertices,basis_type,der)
%����$\int_{x_n}^{x_{n+1}} (u^{der}-w_n)^2 dx,w_n=\sum_{k=1}^{N_{lb}} u_{Tb{k,n}}\fai^{der}_{nk}$
int_value=0;
%Gpn=size(Gauss_point,2);
Gpn=9;
number_of_local_basis = size(Tb,1);
h=1/Nx_max;
for i=1:Gpn
    
    
   
    xable = floor((Gauss_point(1,i)-x_left)/h);
    yable = floor((Gauss_point(2,i)-y_left)/h);
    
    x_data = x_left+(xable)*h;%%�������ұ߽�
    y_data = y_left+(yable)*h;%%
    
    
    if Gauss_point(1,i)+Gauss_point(2,i)<=y_data+x_data+h %%������жϿ����Ǵ���ģ���Ҫ�ú��ٿ�һ��
        position = ((xable)*Ny_max+(yable))*2+1;
        
    else
        position = ((xable)*Ny_max+(yable))*2+2;
        
    end
    
    
            position;
    vertices1=P1(:,T1(:,position));%��ö�Ӧ��Ԫ��ŵĶ�������
    %     [Gauss_coefficient_reference_1D_diff,Gauss_point_reference_1D_diff]=generate_Gauss_reference_triangle(Gauss_type);
    %     [Gauss_weight_diff,Gauss_point_diff]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D_diff,Gauss_point_reference_1D_diff,vertices1);
    %     Gauss_point=Gauss_point_diff';
    u_exact=exact_solution(Tb1(:,position));%��ö�Ӧ��Ԫ�ľֲ���������ϵ�������õ��ģ�
    %����ֲ������ϵ���ֵ���
    
    S = 0;
    for exact_j = 1:number_of_local_basis
        S = S + u_exact(exact_j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices1,basis_type,exact_j,0,0);
    end
    
    
    sum=0;
    for j=1:number_of_local_basis
        %int_value=int_value+Gauss_weight(i)*(feval(exact_solution_fun,Gauss_point(i))-FE_function_1D(Gauss_point(i),uh_local,vertices,basis_type,der))^2;
        sum=sum+uh_local(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,0,0);
    end
    
    int_value=int_value+Gauss_weight(i)*(S-sum)^2;
end
local_error=int_value;
end