function b = assemble_vector_2D_n(solution_n,matrix_size,Pb,Tb,basis_type_trial,basis_type_test)
b=zeros(matrix_size(1),1);
Gpn=9;
number_of_element = size(Tb,2);
number_of_local_basis = size(Tb,1);
for n = 1 : number_of_element
    vertices=Pb(:,Tb(:,n));
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gpn);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    uh_local=solution_n(Tb(:,n));
    for i=1:Gpn
        sum=[];
        sum(i)=0;
        
        for j=1:number_of_local_basis 
            sum(i)=sum(i)+uh_local(j)*FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_trial,j,0,0);
        end
    end
    for beta = 1 : number_of_local_basis
        int_value=0;
        for i=1:Gpn
            int_value=int_value+Gauss_weight(i)*sum(i)* FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type_test,beta,0,0);
        end
        b(Tb(beta,n))=b(Tb(beta,n))+int_value;
    end
end
