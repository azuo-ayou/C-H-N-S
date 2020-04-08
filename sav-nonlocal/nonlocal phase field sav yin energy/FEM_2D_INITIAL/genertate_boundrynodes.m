function boundrynodes = genertate_boundrynodes(Pb,x_left,x_right,y_left,y_right)
%-1:Dirichlet point
%-2:neumann
%-3:robin
% number_of_boundrynodes = 2*(N1+1+N2+1)-4;
% boundrynodes=zeros(2,number_of_boundrynodes);
% boundrynodes(1,:)=-1;
% boundrynodes(2,1)=1;
% boundrynodes(2,2)=size(P,2);
% for k = 1 : number_of_boundrynodes
%     if k <= N1+1
%         boundrynodes(2,k)=(k-1)*(N2+1)+1;
%     else if k <= N1+N2+1
%             boundrynodes(2,k)=N1*(N2+1)+(k-N1);
%         else if k <= (N1+1)*2+N2+1-2
%                 boundrynodes(2,k)=(N1+1-(k-N1-N2-1))*(N2+1);
%             else
%                 boundrynodes(2,k)=N2+1-(k-(2*N1+N2+1));
%             end
%         end
%     end
% end
boundrynodes=[];
number_of_points=size(Pb,2);

for i = 1 : number_of_points
    xx=Pb(1,i);yy=Pb(2,i);
    if xx <= x_left || xx >= x_right || yy <= y_left || yy >= y_right
        boundrynodes=[boundrynodes [-1;i]];
    end
end