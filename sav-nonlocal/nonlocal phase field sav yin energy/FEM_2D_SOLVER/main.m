clear all
clc
format shortE
tic
%    NN = [4,8,16,32];
%  dt = [1/8,1/32,1/64,1/256,1/512];
%  NN = [1/4,1/8,1/16];
%   dt = [1/4,1/8,1/16,1/32];
%  
  NN = [8,16];
 dt = [1/8,1/16];
 
error_u = zeros(3,size(NN,2));
error_p = zeros(3,size(NN,2));
for i=1:size(NN,2)
    format shortE
    %Two_D_FE_solver(202,202,NN(i));
    [error_u(:,i),error_p(:,i)]=Two_D_FE_solver(201,201,NN(i),dt(i));  
end

order_u = zeros(3,size(NN,2));
order_p = zeros(3,size(NN,2));

for i = 1:3
    for j = 1:size(NN,2)-1
        order_u(i,j+1)=log(error_u(i,j)/error_u(i,j+1))/log(2);
        order_p(i,j+1)=log(error_p(i,j)/error_p(i,j+1))/log(2);
    end
end

fprintf('result_fai : \n')
fprintf(' \n')
fprintf('        N          无穷范数          order          L2范数          order            H1范数          order        \n');
for j=1:size(NN,2)
fprintf(' %e ',NN(j));fprintf('   %e   ',error_u(1,j));fprintf('  %e',order_u(1,j));fprintf('   %e',error_u(2,j));fprintf('    %e',order_u(2,j));fprintf('    %e',error_u(3,j));fprintf('    %e\n',order_u(3,j));
end
fprintf(' \n')
fprintf('result_miu : \n')
fprintf(' \n')
fprintf('        N          无穷范数          order         L2范数          order            H1范数          order        \n');
for j=1:size(NN,2)
fprintf(' %e ',NN(j));fprintf('   %e   ',error_p(1,j));fprintf('  %e',order_p(1,j));fprintf('  %e',error_p(2,j));fprintf('    %e',order_p(2,j));fprintf('    %e',error_p(3,j));fprintf('    %e\n',order_p(3,j));
end
toc