function A=assemble_matrix_2D(coe_fun,matrix_size,Pb,Tb,basis_type_trial,der_trial_x,der_trial_y,basis_type_test,der_test_x,der_test_y)
A=sparse(matrix_size(1),matrix_size(2));
Gpn=9;%Gauss points number
number_of_element = size(Tb,2);
number_of_local_basis = size(Tb,1);
for n = 1 : number_of_element
    vertices=Pb(:,Tb(:,n));
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gpn);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    for alpha = 1 : number_of_local_basis
        for beta = 1 : number_of_local_basis 
            int_value =Gauss_vol_int_trial_test(coe_fun,Gauss_weight,Gauss_point,vertices,basis_type_trial,alpha,der_trial_x,der_trial_y,basis_type_test,beta,der_test_x,der_test_y,Gpn);
            A(Tb(beta,n),Tb(alpha,n))= A(Tb(beta,n),Tb(alpha,n))+int_value;           
        end        
    end 
end