% function solution =Two_D_FE_solver(basis_type_trial,basis_type_test)
% clear all clc
% Two_D_FE_solver(201,201)
  function [b_u,b_p]=Two_D_FE_solver_bdf2(basis_type_trial,basis_type_test,N,dt)
% function Two_D_FE_solver(basis_type_trial,basis_type_test,N,dt)
b_u = zeros(3,1);
b_p = zeros(3,1);
N1=N;
N2=N1;
 ef = 0.1;
 gram = 1;
a_dt =dt;

delte = 1/dt;
%  string_2=sprintf('%s%01.2f%s%01.1f%s','ef=',ef,'__gram=',gram,'__');
 string_2=sprintf('%s%d%s%d%s%01.2f%s%01.1f%s','N=',N,'__dt=',delte,'__ef=',ef,'__gram=',gram,'__');
x_left=0;x_right=1;
y_left=0;y_right=1;

xx=0:1/N:1;
yy=0:1/N:1;
[x,y]=meshgrid(xx,yy);

[P,T] = generate_PT_2D(N1,N2,x_left,x_right,y_left,y_right);%the left and right vertices of the interval
% [Pb,Tb]=genertate_PbTb(P,T);
[Pb,Tb]=genertate_PbTb(P,T,basis_type_trial,N1,N2,x_left,x_right,y_left,y_right);
Pb;
 boundrynodes=genertate_boundrynodes(Pb,x_left,x_right,y_left,y_right);

matrix_size=[size(Pb,2),size(Pb,2)];
% matrix_size_2=[size(P,2),size(P,2)]; %%压力部分矩阵

Z1 = zeros(size(Pb,2),size(Pb,2));

solution_n1_3 = exact_u1(Pb(1,:),Pb(2,:),0)';

solution_n1_ex = exact_u1(Pb(1,:),Pb(2,:),dt)';
solution_n2 = exact_u2(Pb(1,:),Pb(2,:),dt)';
solution_n1 = exact_u1(Pb(1,:),Pb(2,:),2*dt)';

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
    E_F_or_rn= assemble_E_2D(solution_n1,solution_n1_ex,Pb,Tb,basis_type_trial,111);%常数 外推;
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
    k = k+1;
end
    
    
    
   A_1= [3*M/(2*dt),A
       gram*ef^2*A,-M]; 
for i=2*dt:dt:1-dt
    
    t = i + a_dt;
%       b11=assemble_vector_2D(@function_f1,matrix_size,Pb,Tb,basis_type_test,i);
      b12=assemble_vector_2D(@function_f1,matrix_size,Pb,Tb,basis_type_test,t);
      b_f = b12;
      % sav
      %第n层的能量
    E_F_or_rn= assemble_E_2D(solution_n1,solution_n1_ex,Pb,Tb,basis_type_trial,333);%常数； % bn = (phi^3-phi)/E_F;
    % bn = (phi^3-phi)/E_F; 
    %第n-1层的能量
    E_F_2 = assemble_E_2D(solution_n1,solution_n1,Pb,Tb,basis_type_trial,111);
    E_F_1 = assemble_E_2D(solution_n1_ex,solution_n1_ex,Pb,Tb,basis_type_trial,333);
    
    
    %A for sav
     A_sav  = assemble_matrix_2D_sav(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0,gram,333);%%   111  和 222都可以
    
      % 含有bn的右端项
      L_bn = assemble_vector_2D_sav_bn(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_test,0,0,gram,333);%向量
      
      
      b_sav = [2*M*solution_n1/dt-M*solution_n1_ex/(2*dt)+b_f
         -(4*E_F_2-E_F_1)*L_bn/3+(4*A_sav*solution_n1/3-A_sav*solution_n1_ex/3)/2];
     b_s = b_sav;
   
     A_s = A_1 + [Z1,Z1
          A_sav/2,Z1];
      
%     [A_s,b_s]=treat_Dirichlet_boundry(A_cn,b,boundrynodes,@exact_u1,@exact_u2,Pb,t);
     solution=A_s\b_s;
     
%     solution_n1_3 = solution_n1_ex;
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
%    if mod((t/dt),1)==0
%        A = figure('Color',[1 1 1]);
%        surf(x,y,reshape(solution_n1,N+1,N+1));
%        view(2)
%        shading interp
%        colormap('jet')
%        axis equal
%        axis square
%        
%      
%         frame = getframe(A);
%         im=frame2im(frame);
% %        string_1=sprintf('%d%s',k,'.jpg');
%        string_1=sprintf('%s%d%s',string_2,k,'.png');
% %        saveas(gcf,string_1,'pdf');
%         print(A,'-dpng', '-r600',string_1); 
% %         imwrite(im,string_1,'jpg');
%    end
    k = k+1;
end
A = figure;
plot(time,E_s);
xlabel('t');
ylabel('Energy');
title('energy');
frame = getframe(A);
im=frame2im(frame);
string_1=sprintf('%s','energy.pdf');
saveas(gcf,string_1);
% imwrite(im,string_1,'jpg');

figure;
plot(time,M_s);
xlabel('t');
ylabel('Mass');
title('mass');

frame = getframe(A);
im=frame2im(frame);
string_1=sprintf('%s','mass.pdf');
saveas(gcf,string_1);
% imwrite(im,string_1,'jpg');
t


end
