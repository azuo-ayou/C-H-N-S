function b = assemble_R_1_2D(solution_n1,Pb,Tb,basis_type)
b=0;
Gpn=9;
number_of_element = size(Tb,2);
number_of_local_basis = size(Tb,1);
for n = 1 : number_of_element
    int_value=0;
    vertices=Pb(:,Tb(:,n));
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gpn);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    uh_local1=solution_n1(Tb(:,n));
%     uh_local2=solution_n1_ex(Tb(:,n));
   for i = 1 : Gpn
%        c2 = 0;
    c1=0;
    for j=1:number_of_local_basis
    c1=c1+uh_local1(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,0,0);
%     c2=c2+uh_local2(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,0,0);
    end
%     c3 = (3*c1-c2)/2;
    %int_value=int_value+Gauss_weight(i)*feval(coe_fun,Gauss_point(1,i),Gauss_point(2,i))*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,basis_index_trial,der_trial_x,der_trial_y)* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,der_test_x,der_test_y);
    int_value = int_value+Gauss_weight(i)*(c1^2-1)^2/4;
   end
   b = b + int_value;
end
b = sqrt(b);
end