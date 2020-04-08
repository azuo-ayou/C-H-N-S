function result = FE_fun_2D(x,y,uh_local,vertices,basis_type,der_x,der_y)
result=0;
% uh_local=solution(Tb(:,n));
for k = 1 : nlb
    result = result + uh_local(k)*FE_local_basis_2D(x,y,vertices,basis_type,k,der_x,der_y);
    
end