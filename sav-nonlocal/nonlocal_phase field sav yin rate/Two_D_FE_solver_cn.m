% function solution =Two_D_FE_solver(basis_type_trial,basis_type_test)
% clear all clc
% Two_D_FE_solver(201,201)
function [b_u,b_p]=Two_D_FE_solver_cn(basis_type_trial,basis_type_test,N,dt)
% function Two_D_FE_solver(basis_type_trial,basis_type_test,N,dt)
b_u = zeros(3,1);
b_p = zeros(3,1);
N1=N;
N2=N1;
% h = 1/N1;
% rou = 1;
%  niu = 1;
ef = 0.1;
gram = 1;
%aph=1;
% time=1;
% dt=16*(1/N)^3;
% dt = 0.2;
% dt=4*(2/N)^2;
% t=1;
a_dt =dt;
x_left=0;x_right=1;
y_left=0;y_right=1;


[P,T] = generate_PT_2D(N1,N2,x_left,x_right,y_left,y_right);%the left and right vertices of the interval
% [Pb,Tb]=genertate_PbTb(P,T);
[Pb,Tb]=genertate_PbTb(P,T,basis_type_trial,N1,N2,x_left,x_right,y_left,y_right);
Pb;
boundrynodes=genertate_boundrynodes(Pb,x_left,x_right,y_left,y_right);

matrix_size=[size(Pb,2),size(Pb,2)];
% matrix_size_2=[size(P,2),size(P,2)]; %%压力部分矩阵

Z1 = zeros(size(Pb,2),size(Pb,2));
% Z2 = zeros(size(P,2),size(P,2));
% Z3 = zeros(size(Pb,2),size(P,2));


%Z2 = zero(size(Pb,2),size(Pb,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution_n1 = zeros(size(Pb,2),1);
% solution_n2 = zeros(size(Pb,2),1);
solution_n1_ex = exact_u1(Pb(1,:),Pb(2,:),0)';
% solution_n2_ex = exact_u2(Pb(1,:),Pb(2,:),0)';
% % solution_n3_ex = exact_p(P(1,:),P(2,:),0)';
solution_n1 = exact_u1(Pb(1,:),Pb(2,:),dt)';
%  solution_n2 = exact_u2(Pb(1,:),Pb(2,:),dt)';
% % solution_n3 = exact_p(P(1,:),P(2,:),0)';
% % size(solution_n)
% solution_L = [solution_n1
%     solution_n2];

%% 矩阵组装 begin

A1=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,1,0,basis_type_test,1,0);
A2=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,0,1,basis_type_test,0,1);
A=(A1+A2);

M=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0);%%质量矩阵


A_1= [M/dt,A
    gram*ef^2*A,-M];
t = 0;
k = 1;
for i=0:dt:0
    t = i + a_dt;
    %     b3 = zeros(size(Pb,2),1);
    %     b_f1 = b3;
    %     b_f2 = b3;
    E_F_or_rn= assemble_E_2D(solution_n1,solution_n1_ex,Pb,Tb,basis_type_trial,111);%常数 外推;type_scheme 111 隐式，222 cn ， 333 bdf2
    E_F= assemble_E_2D(solution_n1,solution_n1,Pb,Tb,basis_type_trial,111);%常数；就是第n项的rn不用外推
    % bn = (phi^3-phi)/E_F;
    
    %A for sav
    A_sav  = assemble_matrix_2D_sav(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0,gram,111);
    
    % 含有bn的右端项
    L_bn = assemble_vector_2D_sav_bn(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_test,0,0,gram,111);%向量
    
    b_sav = [M*solution_n1/dt
        -E_F*L_bn+A_sav*solution_n1/2];
    
    b_s = b_sav;
    
    A_s = A_1 + [Z1,Z1
        A_sav/2,Z1];
    
    % [A_s,b_s]=treat_Dirichlet_boundry(A_cn,b,boundrynodes,@exact_u1,@exact_u2,Pb,t);
    solution=A_s\b_s;
    
    solution_n1_ex = solution_n1;
    
    solution_n1 = solution(1:size(Pb,2));
    solution_n2 = solution(size(Pb,2)+1:size(Pb,2)+size(Pb,2));
    
    E_b3 = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type_trial,ef,gram);
    M_b1 = assemble_M_1_2D(solution_n1,Pb,Tb,basis_type_trial);
    time(k)=t;
    E_s (k)= E_b3;
    M_s (k) =  M_b1;
end






A_1= [M/dt,A/2
    gram*ef^2*A,-M/2];
for i=dt:dt:5-dt
    
    t = i + a_dt;
    
    
    E_F_or_rn= assemble_E_2D(solution_n1,solution_n1_ex,Pb,Tb,basis_type_trial,222);%常数 外推;
    E_F= assemble_E_2D(solution_n1,solution_n1,Pb,Tb,basis_type_trial,222);%常数；就是第n项的rn不用外推
    % bn = (phi^3-phi)/E_F;
    
    %A for sav
    A_sav  = assemble_matrix_2D_sav(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0,gram,111);%% 只能111
    
    % 含有bn的右端项
    L_bn = assemble_vector_2D_sav_bn(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_test,0,0,gram,111);%向量
    
    
%     b_sav = [M*solution_n1/dt-A*solution_n2/2
%         M*solution_n2/2-gram*ef^2*A*solution_n1/2-E_F*L_bn+A_sav*solution_n1/4];
%     
    b_sav = [M*solution_n1/dt-A*solution_n2/2
        M*solution_n2/2-E_F*L_bn+A_sav*solution_n1/4];
  
    
    b_s = b_sav;
    
    A_s = A_1 + [Z1,Z1
        A_sav/4,Z1];
    
    %     [A_s,b_s]=treat_Dirichlet_boundry(A_cn,b,boundrynodes,@exact_u1,@exact_u2,Pb,t);
    solution=A_s\b_s;
    
    solution_n1_ex = solution_n1;
    
    solution_n1 = solution(1:size(Pb,2));
    solution_n2 = solution(size(Pb,2)+1:size(Pb,2)+size(Pb,2));
    
    E_b3 = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type_trial,ef,gram);
    M_b1 = assemble_M_1_2D(solution_n1,Pb,Tb,basis_type_trial);
    %      R_b2 = assemble_R_1_2D(solution_n1,Pb,Tb,basis_type_trial);
    time(k)=t;
    %     E_s (k)= (ef^2*E_b1+E_b2)/(4*dt)+(R_b2+(2*R_b2-R_b1)^2)/(2*dt);
    E_s (k)= E_b3;
    M_s (k) =  M_b1;
    k = k+1;
end

 A = figure;
plot(time,E_s);
xlabel('t');
ylabel('Enery');
% title('energy');
% frame = getframe(A);
% im=frame2im(frame);
% string_1=sprintf('%s','energy.jpg');
% imwrite(im,string_1,'jpg');

 figure;
plot(time,M_s);
xlabel('t');
ylabel('Mass');
% title('mass');
% frame = getframe(A);
% im=frame2im(frame);
% string_1=sprintf('%s','mass.jpg');
% imwrite(im,string_1,'jpg');
t

end
