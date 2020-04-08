format shortE
clear all
% 计算收敛阶请先画图比较
basis_type = 201;

% NN = [30,30];
% dt = [1/30,1/30];

x_left=-0.5;x_right=0.5;
y_left=-0.5;y_right=0.5;

NN = [10,15,20,25,30,35,40,45,50,100];
dt = [1/10,1/15,1/20,1/25,1/30,1/35,1/40,1/45,1/50,1/100];

% NN = [2,4,8,16,32,64];
% dt = [1/2,1/4,1/8,1/16,1/32,1/64];

% error_u = zeros(3,size(NN,2));
% error_p = zeros(3,size(NN,2));
S = size(NN,2);
[solution2,Tb,T,P]= Two_D_FE_solver_yin(basis_type,basis_type,NN(S),dt(S));
%  k=1;
for i=1:S-1
    format shortE
    %      solution1 = Two_D_FE_solver(basis_type,basis_type,NN(i),dt(i));
    [solution1,Tb1,T1,P1]= Two_D_FE_solver_yin(basis_type,basis_type,NN(i),dt(i));
    
%     size(error)
    
%     err_L2(i) = L2_norm_2D(error,@exact_0,Tb1,T1,P1,basis_type,1);
    err_L2(i) = L2_norm_2D_interpolation_error_int(solution1,solution2,Tb1,T1,P1,basis_type,1,NN(S),NN(S),Tb,T,P,x_left,y_left);

end
%                           err;
for j = 1:size(err_L2,2)-1
    rate(j) = (log(err_L2(j)/err_L2(j+1)))/(log(NN(j+1)/NN(j)));
end

loglog(dt(1:end-2),err_L2(1:end-1),'-ro','Linewidth',2);
% hold on
% loglog(dt,error_u(2,:),'-Ko','Linewidth',2);
% hold on
% loglog(dt,error_u(3,:),'-Go','Linewidth',2);

