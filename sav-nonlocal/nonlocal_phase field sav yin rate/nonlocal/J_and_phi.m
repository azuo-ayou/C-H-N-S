function A=J_and_phi(Exact_J,matrix_size,Pb,Tb,basis_type_trial,basis_type_test)
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
            int_value =Gauss_vol_int_trial_test_AN_nonlocal(Exact_J,Gauss_weight,Gauss_point,vertices,basis_type_trial,alpha,basis_type_test,beta,Gpn);
            A(Tb(beta,n),Tb(alpha,n))= A(Tb(beta,n),Tb(alpha,n))+int_value;           
        end        
    end 
end

function int_value = Gauss_vol_int_trial_test_AN_nonlocal(Exact_J,Gauss_weight,Gauss_point,vertices,basis_type_trial,basis_index_trial,basis_type_test,basis_index_test,Gpn)
int_value=0;
for k= 1 : Gpn
    int_value1=0;
    for i = 1 : Gpn
%         int_value=int_value+Gauss_weight(i)*feval(Exact_J,Gauss_point(1,k)-Gauss_point(1,i),Gauss_point(2,k)-Gauss_point(2,i))*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,basis_index_trial,0,0)* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,0,0);

        int_value1=int_value1+Gauss_weight(i)*feval(Exact_J,Gauss_point(1,k)-Gauss_point(1,i),Gauss_point(2,k)-Gauss_point(2,i))*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,basis_index_trial,0,0);
    end
    int_value=int_value+int_value1* Gauss_weight(k)*FE_local_basis_2D(Gauss_point(1,k),Gauss_point(2,k),vertices,basis_type_test,basis_index_test,0,0);
end