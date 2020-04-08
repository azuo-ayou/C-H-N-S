N = 16;
xx=0:1/N:1;
yy=0:1/N:1;
[x,y]=meshgrid(xx,yy);
figure;
% surf(x,y,reshape(solution_n1,N+1,N+1));
surf(x,y,reshape(exact_u1(Pb(1,:),Pb(2,:),0)',N+1,N+1));