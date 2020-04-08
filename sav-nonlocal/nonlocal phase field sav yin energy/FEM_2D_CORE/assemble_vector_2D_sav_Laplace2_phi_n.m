function b = assemble_vector_2D_sav_Laplace2_phi_n(fb1,matrix_size,Pb,Tb,basis_type_test)
b=zeros(matrix_size(1),1);
Gpn=9;
number_of_element = size(Tb,2);
number_of_local_basis = size(Tb,1);
for n = 1 : number_of_element
    vertices=Pb(:,Tb(:,n));
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gpn);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    
    uh_local_1=fb1(Tb(:,n));
    
    for beta = 1 : number_of_local_basis
        int_value = Gauss_vol_int_f_bN(uh_local_1,number_of_local_basis,Gauss_weight,Gauss_point,vertices,basis_type_test,beta,Gpn);
        b(Tb(beta,n))=b(Tb(beta,n))+int_value;
    end
    
end

function int_value = Gauss_vol_int_f_bN(uh_local_1,number_of_local_basis,Gauss_weight,Gauss_point,vertices,basis_type_test,basis_index_test,Gpn)
int_value=0;
for i = 1 : Gpn
    c1=0;
    c2=0;
%     c3=0;
    for j=1:number_of_local_basis
    %直接求拉普拉斯
    c1=c1+uh_local_1(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,j,2,0);
    c2=c2+uh_local_1(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,j,0,2);

    
    %转移两次导数
%     c1=c1+uh_local_1(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,j,AN_der1_x,AN_der1_y);
    end
    
      int_value=int_value+Gauss_weight(i)*(c1+c2)*(FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,2,0)+FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,basis_index_test,0,2));

end