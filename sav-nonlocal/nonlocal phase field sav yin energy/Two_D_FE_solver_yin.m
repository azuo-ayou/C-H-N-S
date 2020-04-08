% function solution =Two_D_FE_solver(basis_type_trial,basis_type_test)
% clear all clc
% Two_D_FE_solver(201,201)
function [b_u,b_p]=Two_D_FE_solver_yin(basis_type_trial,basis_type_test,N,dt)
% function Two_D_FE_solver(basis_type_trial,basis_type_test,N,dt)
b_u = zeros(3,1);
b_p = zeros(3,1);
N1=N;
N2=N1;
% ef = 0.01;
% ef_1 = 1;
ef = 0.02;
ef_1 = 1;
gram =1;
beta = 2;
delte = 1/dt;
string_2=sprintf('%s%d%s%d%s%01.2f%s%01.1f%s','N=',N,'__dt=',delte,'__ef=',ef,'__gram=',gram,'__');

period = 8;
Time = 2;
a_dt =dt;
t=0;
k=1;
x_left=-0.5;x_right=0.5;
y_left=-0.5;y_right=0.5;
[P,T] = generate_PT_2D(N1,N2,x_left,x_right,y_left,y_right);%the left and right vertices of the interval
[Pb,Tb]=genertate_PbTb(P,T,basis_type_trial,N1,N2,x_left,x_right,y_left,y_right);
Pb;
% boundrynodes=genertate_boundrynodes(Pb,x_left,x_right,y_left,y_right);
matrix_size=[size(Pb,2),size(Pb,2)];
Z1 = zeros(size(Pb,2),size(Pb,2));
solution_n1_ex = exact_u1(Pb(1,:),Pb(2,:),0)';
solution_n1 = exact_u1(Pb(1,:),Pb(2,:),dt)';
solution_n2 = exact_u2(Pb(1,:),Pb(2,:),dt)';

xx=-0.5:1/N:0.5;
yy=-0.5:1/N:0.5;
[x,y]=meshgrid(xx,yy);

E_b3 = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type_trial,ef,gram,beta);
M_b1 = assemble_M_1_2D(solution_n1,Pb,Tb,basis_type_trial);
time(k)=t;
E_s (k)= E_b3;
M_s (k) =  M_b1;
Ap = figure('Color',[1 1 1]);
surf(x,y,reshape(solution_n1,N+1,N+1));
view(2)
shading interp
colormap('jet')
axis square
frame = getframe(Ap);
im=frame2im(frame);

string_1=sprintf('%s%d%s',string_2,k,'.jpg');
imwrite(im,string_1,'jpg');
k = k+1;




A1=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,1,0,basis_type_test,1,0);
A2=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,0,1,basis_type_test,0,1);
A=(A1+A2);
M=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0);%%质量矩阵


J0=int_fun(@Exact_J,Pb,Tb);%计算J的积分
J_phi_1=J_and_phi(@Exact_J,matrix_size,Pb,Tb,basis_type_trial,basis_type_test);%计算J和phi的两重积分
J_phi_2=J_and_phi_add_xx(@Exact_J,matrix_size,Pb,Tb,basis_type_trial,basis_type_test);
J_phi=J_phi_2-J_phi_1;
J_SUMMARY = J_and_phi_only_once(@Exact_J,matrix_size,Pb,Tb,basis_type_trial,basis_type_test);

max(max(J_SUMMARY-J_phi))
min(min(J_SUMMARY-J_phi))
J_phi = J_SUMMARY;






A_1= [M/dt,J_phi
    gram*ef^2*J_phi+beta*gram*M/ef_1^2,-M];

% A_1= [M/dt,A
%     gram*ef^2*A+beta*gram*M/ef_1^2,-M];


type_scheme = 111;
for i=0:dt:Time-dt
    t = i + a_dt;
    
    E_F_or_rn= assemble_E_2D(solution_n1,solution_n1_ex,Pb,Tb,basis_type_trial,type_scheme,beta);%常数 外推;
    E_F= assemble_E_2D(solution_n1,solution_n1,Pb,Tb,basis_type_trial,type_scheme,beta);%常数；就是第n项的rn不用外推
    
    
    %A for sav
    A_sav  = assemble_matrix_2D_sav(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0,gram,type_scheme,beta);
    
    % 含有bn的右端项
    L_bn = assemble_vector_2D_sav_bn(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_test,0,0,gram,type_scheme,beta);%向量
    
    b_sav = [M*solution_n1/dt
        (-E_F*L_bn+A_sav*solution_n1/2)/ef_1^2];
    
    b_s = b_sav;
    
    A_s = A_1 + [Z1,Z1
        A_sav/2/ef_1^2,Z1];
    
    % [A_s,b_s]=treat_Dirichlet_boundry(A_cn,b,boundrynodes,@exact_u1,@exact_u2,Pb,t);
    solution=A_s\b_s;
    
    solution_n1_ex = solution_n1;
    
    solution_n1 = solution(1:size(Pb,2));
    solution_n2 = solution(size(Pb,2)+1:size(Pb,2)+size(Pb,2));
    
    E_b3 = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type_trial,ef,gram,beta);
    M_b1 = assemble_M_1_2D(solution_n1,Pb,Tb,basis_type_trial);
    time(k)=t;
    E_s (k)= E_b3;
    M_s (k) =  M_b1;
    
    
    if mod((t/dt),period)==0
        Ap = figure('Color',[1 1 1]);
        surf(x,y,reshape(solution_n1,N+1,N+1));
        view(2)
        shading interp
        colormap('jet')
        axis equal
        axis square
        
        frame = getframe(Ap);
        im=frame2im(frame);
        %        string_1=sprintf('%d%s',k,'.jpg');
        string_1=sprintf('%s%d%s',string_2,k,'.jpg');
        imwrite(im,string_1,'jpg');
    end
    
    k = k+1;
    
    
end

Ap = figure;
plot(time,E_s);
xlabel('t');
ylabel('Energy');
title('energy');
frame = getframe(Ap);
im=frame2im(frame);
string_1=sprintf('%s','energy.jpg');
imwrite(im,string_1,'jpg');

Ap = figure;
plot(time,M_s);
xlabel('t');
ylabel('Mass');
title('mass');
frame = getframe(Ap);
im=frame2im(frame);
string_1=sprintf('%s','mass.jpg');
imwrite(im,string_1,'jpg');
t


end
