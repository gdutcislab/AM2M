function p=MaOP12(p)
 p.name='MaOP12';
 p.od = 10;
 p.pd = 10;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain        =  range;
 p.func          = @evaluate;
    %P13  evaluation function.
    function ff = evaluate(x)
        [dim, num]   = size(x);
        g            = sum((x(3:dim,:) - (x(1,:)+x(2,:))/2).^2);       
        h            = zeros(size(x,2),1);
        for i=1:num
           if abs(x(1,i)*x(2,i)- x(1,i)*(1- x(2,i)))>0.001
              h(i)   =  max(0,-(1+0.4)*sin(2*pi*(x(1,i))));           
           end
        end
         y(1,:)       = x(1,:).*x(2,:);
         y(2,:)       = x(1,:).*(1-x(2,:));
         y(3,:)       = (1-x(1,:));
         ff(1:3,:)    = (1+g+h).*y(1:3,:);
         ff(4,:)      = (0.1*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(4,:)-(x(1,:)+x(2,:))/2).^2)-1);
         ff(5,:)      = (0.4*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(5,:)-(x(1,:)+x(2,:))/2).^2)-1);
         ff(6,:)      = (0.3*y(1,:)+0.2*y(2,:)+0.4*y(3,:)).*(exp((x(6,:)-(x(1,:)+x(2,:))/2).^2)-1);
         ff(7,:)      = (0.3*y(1,:)+0.1*y(2,:)+0.1*y(3,:)).*(exp((x(7,:)-(x(1,:)+x(2,:))/2).^2)-1);
         ff(8,:)      = (0.1*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(8,:)-(x(1,:)+x(2,:))/2).^2)-1);
         ff(9,:)      = (0.1*y(1,:)+0.2*y(2,:)+0.3*y(3,:)).*(exp((x(9,:)-(x(1,:)+x(2,:))/2).^2)-1);
         ff(10,:)     = (0.1*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(10,:)-(x(1,:)+x(2,:))/2).^2)-1); 
    end
end