function result = function_g1(x,y,t)
%the Dirichlet boundry function 
 x_left=0;x_right=1;
 y_left=-0.25;y_right=0;
% result = 0;
if x == x_left
    result =(exp(-y)).*cos(2*pi*t);
else if x == x_right
        result = (y.^2+exp(-y))*cos(2*pi*t);
    else if y == y_left
            result = (1/16*x.^2+exp(0.25)).*cos(2*pi*t);
        else if y == y_right
                result = cos(2*pi*t) ;
%             else 
%                 warning('input (x,y) is not a Dirichlet point');
             end
        end
    end
end
