function p=MaOP7(p)
 p.name='MaOP7';
 p.od = 10;
 p.pd = 10;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %DTLZ1  evaluation function.
    function ff = evaluate(x)
          [dim, num]   = size(x);
          y            = zeros(p.od,num);
          y(1,:)       = x(1,:);
          y(2,:)       = x(2,:);
          temp         = sum((x(3:dim,:)-x(1,:).*x(2,:)).^2);
          y(3,:)       = 2-sum(y(1:2,:).*(1+sin(3*pi*y(1:2,:)))/2);        
          ff(1:3,:)    = (1+temp).*y(1:3,:);
          ff(4,:)      = (0.1*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(4,:)-x(1,:).*x(2,:)).^2)-1);
          ff(5,:)      = (0.4*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(5,:)-x(1,:).*x(2,:)).^2)-1);
          ff(6,:)      = (0.3*y(1,:)+0.2*y(2,:)+0.4*y(3,:)).*(exp((x(6,:)-x(1,:).*x(2,:)).^2)-1);
          ff(7,:)      = (0.3*y(1,:)+0.1*y(2,:)+0.1*y(3,:)).*(exp((x(7,:)-x(1,:).*x(2,:)).^2)-1);
          ff(8,:)      = (0.1*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(8,:)-x(1,:).*x(2,:)).^2)-1);
          ff(9,:)      = (0.1*y(1,:)+0.2*y(2,:)+0.3*y(3,:)).*(exp((x(9,:)-x(1,:).*x(2,:)).^2)-1);
          ff(10,:)     = (0.1*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(10,:)-x(1,:).*x(2,:)).^2)-1); 
%          ff(4,:)       = (0.1*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(4,:)-0.5).^2)-1);
%          ff(5,:)       = (0.4*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(5,:)-0.5).^2)-1);
%          ff(6,:)       = (0.3*y(1,:)+0.2*y(2,:)+0.4*y(3,:)).*(exp((x(6,:)-0.5).^2)-1);
%          ff(7,:)       = (0.3*y(1,:)+0.1*y(2,:)+0.1*y(3,:)).*(exp((x(7,:)-0.5).^2)-1);
%          ff(8,:)       = (0.1*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(8,:)-0.5).^2)-1);
%          ff(9,:)       = (0.1*y(1,:)+0.2*y(2,:)+0.3*y(3,:)).*(exp((x(9,:)-0.5).^2)-1);
%          ff(10,:)      = (0.1*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(10,:)-0.5).^2)-1);
    end
end