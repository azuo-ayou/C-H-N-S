function r = int_fun(func,Pb,Tb)
% A=sparse(matrix_size(1),matrix_size(2));
Gpn=9;%Gauss points number
number_of_element = size(Tb,2);
r =0;
% number_of_local_basis = size(Tb,1);
for n = 1 : number_of_element
    vertices=Pb(:,Tb(:,n));
    [Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_triangle(Gpn);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices);
    Gauss_point=Gauss_point';
    int_value1 = 0;
    for i = 1 : Gpn
        int_value1=int_value1+Gauss_weight(i)*feval(func,Gauss_point(1,i),Gauss_point(2,i));
    end
    r = r + int_value1;
end