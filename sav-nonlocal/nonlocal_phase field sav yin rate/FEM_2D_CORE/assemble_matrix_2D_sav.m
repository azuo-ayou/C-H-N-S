function A=assemble_matrix_2D_sav(a,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_trial,der_trial_x,der_trial_y,basis_type_test,der_test_x,der_test_y,gram,type_scheme,beta)
A=sparse(matrix_size(1),matrix_size(2));
Gpn=9;%Gauss points number
number_of_element = size(Tb,2);
number_of_local_basis = size(Tb,1);
for n = 1 : number_of_element
    vertices=Pb(:,Tb(:,n));
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gpn);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    uh_local1=solution_n1(Tb(:,n));
    uh_local2=solution_n1_ex(Tb(:,n));
    for alpha = 1 : number_of_local_basis
        for beta = 1 : number_of_local_basis 
            int_value =Gauss_vol_int_trial_test_sav(number_of_local_basis,a,uh_local1,uh_local2,Gauss_weight,Gauss_point,vertices,basis_type_trial,alpha,der_trial_x,der_trial_y,basis_type_test,beta,der_test_x,der_test_y,Gpn,gram,type_scheme,beta);
            A(Tb(beta,n),Tb(alpha,n))= A(Tb(beta,n),Tb(alpha,n))+int_value;           
        end        
    end 
end
end

function int_value = Gauss_vol_int_trial_test_sav(number_of_local_basis,a,uh_local1,uh_local2,Gauss_weight,Gauss_point,vertices,basis_type_trial,basis_index_trial,der_trial_x,der_trial_y,basis_type_test,basis_index_test,der_test_x,der_test_y,Gpn,gram,type_scheme,beta)
int_value1=0;
int_value2=0;
for i = 1 : Gpn
    c1=0;
    c2 =0;
    for j=1:number_of_local_basis
    c1=c1+uh_local1(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,j,0,0);
    c2=c2+uh_local2(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,j,0,0);
    end
    
   if type_scheme == 111
       c3 = c1;
   elseif type_scheme == 222
      c3 = (3*c1-c2)/2;
   elseif type_scheme == 333
      c3 = (2*c1-c2);
   end
    %int_value=int_value+Gauss_weight(i)*feval(coe_fun,Gauss_point(1,i),Gauss_point(2,i))*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,basis_index_trial,der_trial_x,der_trial_y)* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,der_test_x,der_test_y);

int_value1=int_value1+Gauss_weight(i)*gram*(c3^3-c3-beta*c3)/a* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,basis_index_trial,der_trial_x,der_trial_y);
int_value2=int_value2+Gauss_weight(i)*gram*(c3^3-c3-beta*c3)/a* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,der_test_x,der_test_y);

%     int_value=int_value+Gauss_weight(i)*feval(coe_fun,Gauss_point(1,i),Gauss_point(2,i))*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,basis_index_trial,der_trial_x,der_trial_y)* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,der_test_x,der_test_y);
                         
end
int_value=int_value1*int_value2;
end