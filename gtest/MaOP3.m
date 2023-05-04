function p=MaOP3(p)
 p.name='MaOP3';
 p.od = 10;
 p.pd = 10;%这里可以改为4，如果要求自变量空间是4维的情况
 range      = ones(p.pd,2); 
 range(:,1) = 0; 
 p.domain   =  range;
 p.func     = @evaluate;
    %MaOP3  evaluation function.
     function ff = evaluate(x)  
         [dim, num]   = size(x);    
         g            = sum((x(2:dim,:) - sin(pi*repmat(x(1,:),[dim-1,1])/2)).^2); 
         y            = zeros(p.od,num);   
         tmp          = max(0,2*sin(4*pi*x(1,:)));
         y(1,:)       = sqrt(2)/2*x(1,:);
         y(2,:)       = sqrt(2)/2*x(1,:); 
         y(3,:)       = 3*(1.0 - x(1,:).^2); 
         ff(1:3,:)    = (1+tmp+g).*y(1:3,:);
         ff(4,:)      = (0.1*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(4,:)-sin(pi*x(1,:)/2)).^2)-1);
         ff(5,:)      = (0.4*y(1,:)+0.1*y(2,:)+0.3*y(3,:)).*(exp((x(5,:)-sin(pi*x(1,:)/2)).^2)-1);
         ff(6,:)      = (0.3*y(1,:)+0.2*y(2,:)+0.4*y(3,:)).*(exp((x(6,:)-sin(pi*x(1,:)/2)).^2)-1);
         ff(7,:)      = (0.3*y(1,:)+0.1*y(2,:)+0.1*y(3,:)).*(exp((x(7,:)-sin(pi*x(1,:)/2)).^2)-1);
         ff(8,:)      = (0.4*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(8,:)-sin(pi*x(1,:)/2)).^2)-1);
         ff(9,:)      = (0.3*y(1,:)+0.2*y(2,:)+0.3*y(3,:)).*(exp((x(9,:)-sin(pi*x(1,:)/2)).^2)-1);
         ff(10,:)     = (0.1*y(1,:)+0.3*y(2,:)+0.1*y(3,:)).*(exp((x(10,:)-sin(pi*x(1,:)/2)).^2)-1);       
        clear Y tmp1 tmp2;
     end
end