function result = function_g2(x,y,t)
%the Dirichlet boundry function 
 x_left=0;x_right=1;
 y_left=-0.25;y_right=0;
% result = 0;
if x == x_left
    result =2*cos(2*pi*t);
else if x == x_right
        result = (-2/3*y.^3+2)*cos(2*pi*t);
    else if y == y_left
            result = (1/96*x+2-pi*sin(pi*x))*cos(2*pi*t);
        else if y == y_right
                result = (2-pi*sin(pi*x))*cos(2*pi*t) ;
%             else 
%                 warning('input (x,y) is not a Dirichlet point');
             end
        end
    end
end
