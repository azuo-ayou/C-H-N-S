function b = assemble_bn_with_f_n_2D(a,solution_n1,Pb,Tb,basis_type,f,t)
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
    uh_local=solution_n1(Tb(:,n));
   for i = 1 : Gpn    
    c1=0;
    for j=1:number_of_local_basis
    c1=c1+uh_local(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,0,0);
    end
    %int_value=int_value+Gauss_weight(i)*feval(coe_fun,Gauss_point(1,i),Gauss_point(2,i))*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,basis_index_trial,der_trial_x,der_trial_y)* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,der_test_x,der_test_y);
    int_value = int_value+Gauss_weight(i)*(c1^3-c1)/a*feval(f,Gauss_point(1,i),Gauss_point(2,i),t);
   end
   b = b + int_value;
end
end