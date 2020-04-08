function b = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type,ef,gram,beta)% ||grid(phi_n+1)||^2
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
    
    for k = 1 : Gpn
        c1 = 0;
%         c3 = 0;
        int_value1 = 0;
        
        for i = 1:Gpn
            
            c3 = 0;
            
            for j=1:number_of_local_basis
%                 c3=c3+uh_local1(j)*FE_local_basis_2D(Gauss_point(1,k)-Gauss_point(1,i),Gauss_point(2,k)-Gauss_point(2,i),vertices,basis_type,j,0,0);
                c3=c3+uh_local1(j)*(FE_local_basis_2D(Gauss_point(1,k),Gauss_point(2,k),vertices,basis_type,j,0,0)-FE_local_basis_2D(Gauss_point(1,i),Gauss_point(2,i),vertices,basis_type,j,0,0));
            end
            
            %c2 = feval(Exact_J,(Gauss_point(1,k)-Gauss_point(1,i)),(Gauss_point(2,k)-Gauss_point(2,i)));
            c2=Exact_J(Gauss_point(1,k)-Gauss_point(1,i),(Gauss_point(2,k)-Gauss_point(2,i)));
            int_value1 = int_value1+Gauss_weight(i)*c3^2*c2;
        end
        
        
        for j=1:number_of_local_basis
            c1=c1+uh_local1(j)*FE_local_basis_2D(Gauss_point(1,k),Gauss_point(2,k),vertices,basis_type,j,0,0);      
        end
        
        
        int_value = int_value+Gauss_weight(k)*(gram*ef^2*int_value1/4+gram*(c1^2-1-beta)^2/4+beta*c1^2/2-(beta^2+2*beta)/4);
    end
    
    b = b + int_value;
end
% b = sqrt(b);
end