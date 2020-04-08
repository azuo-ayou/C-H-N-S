% function solution =Two_D_FE_solver(basis_type_trial,basis_type_test)
% clear all clc
% Two_D_FE_solver(201,201)
%   function [b_u,b_p]=Two_D_FE_solver(basis_type_trial,basis_type_test,N,dt)
 function N = Two_D_FE_solver(basis_type_trial,basis_type_test,N,dt,type_scheme)
b_u = zeros(3,1);
b_p = zeros(3,1);
N1=N;
N2=N1;
% h = 1/N1;
% rou = 1;
%  niu = 1;
 ef = 0.07;
 ef_2 = 1;
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
for j = 1:size(Pb,2)
solution_n1_3(j) = exact_u1(Pb(1,j),Pb(2,j),0);
end
solution_n1_3 = solution_n1_3';
exact_u1_0 = solution_n1_3;
solution_n1 = solution_n1_3;
solution_n1_ex = solution_n1_3;

%% 矩阵组装 begin
xx=0:1/N:1;
yy=0:1/N:1;
[x,y]=meshgrid(xx,yy);

    A1=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,1,0,basis_type_test,1,0);
    A2=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,0,1,basis_type_test,0,1);
    A=(A1+A2);
    
    M=assemble_matrix_2D(@function_c,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0);%%质量矩阵

    A_1= [3*M/(2*dt),A
       gram*ef^2*A,-M]; %隐式
%     A_1= [3*M/(2*dt),A
%        A,-M]; %隐式

  A_2 = [M/dt,A
       gram*ef^2*A,-M]; %隐式
    t = 0;
    k = 1;
    
    E_b3 = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type_trial,ef,gram);
M_b1 = assemble_M_1_2D(solution_n1,Pb,Tb,basis_type_trial);
   E_s (k)= E_b3;
   M_s (k) =  M_b1;
   time(k)=0;
    k = k+1; 
A = figure;
 surf(x,y,reshape(solution_n1,N+1,N+1));
 view(2)
 shading interp
 colormap('jet')
 axis equal
 frame = getframe(A);
im=frame2im(frame);
string_1=sprintf('%d%s',k,'.jpg');
imwrite(im,string_1,'jpg');
%     b_f=assemble_vector_2D(@function_f1,matrix_size,Pb,Tb,basis_type_test,dt);
     b_f = 0;
%     E_F_or_rn= assemble_E_2D(solution_n1,solution_n1_ex,Pb,Tb,basis_type_trial);%常数； % bn = (phi^3-phi)/E_F;
    % bn = (phi^3-phi)/E_F;  
   %第n-1层的能量
    E_F = assemble_E_2D(solution_n1_ex,solution_n1_3,Pb,Tb,basis_type_trial,type_scheme);   
    %A for sav
     A_sav  = assemble_matrix_2D_sav(E_F,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0,gram,type_scheme);
      % 含有bn的右端项
      L_bn = assemble_vector_2D_sav_bn(E_F,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_test,0,0,gram,type_scheme);%向量
      b_sav = [M*solution_n1/dt+b_f
         (-(E_F)*L_bn+A_sav*solution_n1/2)/ef_2^2];
     b = b_sav;
   
     A_cn = A_2 + [Z1,Z1
          A_sav/2,Z1];
      
%       [A_cn,b]=treat_Dirichlet_boundry(A_cn,b,boundrynodes,@exact_u1,@exact_u2,Pb,t);
     solution=A_cn\b;
     
    solution_n1_3 = solution(1:size(Pb,2));
    solution_n1_ex = solution(1:size(Pb,2));
     
    solution_n1 = solution(1:size(Pb,2));
%     solution_n2 = solution(size(Pb,2)+1:size(Pb,2)+size(Pb,2));
    
    E_b3 = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type_trial,ef,gram);
   M_b1 = assemble_M_1_2D(solution_n1,Pb,Tb,basis_type_trial);
   E_s (k)= E_b3;
   M_s (k) =  M_b1;
   time(k)=dt;
    k = k+1; 
 A=    figure;
 surf(x,y,reshape(solution_n1,N+1,N+1));
%    surf(x,y,reshape(solution_n1,N+1,N+1));
 view(2)
 shading interp
 colormap('jet')
  frame = getframe(A);
im=frame2im(frame);
string_1=sprintf('%d%s',k,'.jpg');
imwrite(im,string_1,'jpg');  

   
for i=dt:dt:1-dt
    
    t = i + a_dt;
      b_f = 0;
    E_F_or_rn= assemble_E_2D(solution_n1,solution_n1_ex,Pb,Tb,basis_type_trial,type_scheme);%常数； % bn = (phi^3-phi)/E_F;
    % bn = (phi^3-phi)/E_F;
    
    %第n-1层的能量
    E_F = assemble_E_2D(solution_n1,solution_n1,Pb,Tb,basis_type_trial,type_scheme);%n
    E_F_n = assemble_E_2D(solution_n1_ex,solution_n1_ex,Pb,Tb,basis_type_trial,type_scheme);%n-1
  
    %A for sav
     A_sav  = assemble_matrix_2D_sav(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_trial,0,0,basis_type_test,0,0,gram,type_scheme);
    
      % 含有bn的右端项
      L_bn = assemble_vector_2D_sav_bn(E_F_or_rn,solution_n1,solution_n1_ex,matrix_size,Pb,Tb,basis_type_test,0,0,gram,type_scheme);%向量
      
      
      b_sav = [2*M*solution_n1/dt-M*solution_n1_ex/(2*dt)+b_f
         (-(4*E_F-E_F_n)*L_bn/3+4*A_sav*solution_n1/3-A_sav*solution_n1_ex/3)/ef_2^2];
     b = b_sav;
   
     A_cn = A_1 + [Z1,Z1
          A_sav,Z1];
      
%       [A_cn,b]=treat_Dirichlet_boundry(A_cn,b,boundrynodes,@exact_u1,@exact_u2,Pb,t);
     solution=A_cn\b;
%       R_b1 = assemble_R_1_2D(solution_n1,Pb,Tb,basis_type_trial);
    solution_n1_3 = solution_n1_ex;
    solution_n1_ex = solution_n1;
     
    solution_n1 = solution(1:size(Pb,2));
%     solution_n2 = solution(size(Pb,2)+1:size(Pb,2)+size(Pb,2));
   
      E_b3 = assemble_E_3_2D(solution_n1,Pb,Tb,basis_type_trial,ef,gram);
     M_b1 = assemble_M_1_2D(solution_n1,Pb,Tb,basis_type_trial);
%      R_b2 = assemble_R_1_2D(solution_n1,Pb,Tb,basis_type_trial);
    time(k)=t;
%     E_s (k)= (ef^2*E_b1+E_b2)/(4*dt)+(R_b2+(2*R_b2-R_b1)^2)/(2*dt);
    E_s (k)= E_b3;
   M_s (k) =  M_b1;
    k = k+1;
%     if mod((t/dt),10)==0
%   A =  figure;
%  surf(x,y,reshape(solution_n1,N+1,N+1));
%  view(2)
%  shading interp
%  colormap('jet')
%  frame = getframe(A);
% im=frame2im(frame);
% string_1=sprintf('%d%s',k,'.jpg');
% imwrite(im,string_1,'jpg');
%     end
%  title('')
end
A = figure;
plot(time,E_s);
title('energy');
frame = getframe(A);
im=frame2im(frame);
string_1=sprintf('%s','energy.jpg');
imwrite(im,string_1,'jpg');

figure;
plot(time,M_s);
title('mass');
frame = getframe(A);
im=frame2im(frame);
string_1=sprintf('%s','mass.jpg');
imwrite(im,string_1,'jpg');
    
    
    figure;
 surf(x,y,reshape(solution_n1,N+1,N+1));   
 view(2)
 shading interp
 colormap('jet')

t

 end
