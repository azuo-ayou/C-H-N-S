function b = assemble_vector_2D(fb,matrix_size,Pb,Tb,basis_type_test,t)
b=zeros(matrix_size(1),1);
Gpn=9;
number_of_element = size(Tb,2);
number_of_local_basis = size(Tb,1);
for n = 1 : number_of_element
    vertices=Pb(:,Tb(:,n));
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gpn);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    for beta = 1 : number_of_local_basis
        int_value = Gauss_vol_int_f(fb,Gauss_weight,Gauss_point,t,vertices,basis_type_test,beta,Gpn);
        b(Tb(beta,n))=b(Tb(beta,n))+int_value;
    end
end
