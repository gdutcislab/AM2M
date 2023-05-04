function pf = pareto(name, no, dim)

    if nargin<3, dim = 3; end
    if nargin<2, no  = 300; end
    switch name
        case 'MaOP1'           
            num           = floor(no/2);
            pf            = zeros(3,no);            
            pf(1,1:num)   = sqrt(2)/2*linspace(0,0.25,num);
            pf(1,num+1:no)= sqrt(2)/2*linspace(0.75,1.0,no-num); 
            pf(2,:)       = pf(1,:);
            pf(3,:)       = 1-(sqrt(2)*pf(1,:)).^2;                       
        case 'MaOP2'
            num           = floor(no/2);
            pf            = zeros(3,no);            
            pf(1,1:num)   = sqrt(2)/2*linspace(0,0.25,num);
            pf(1,num+1:no)= sqrt(2)/2*linspace(0.75,1.0,no-num); 
            pf(2,:)       = pf(1,:);            
            pf(3,:)       = 1-(sqrt(2)*pf(1,:)).^(1/2);         
        case 'MaOP3'
            num           = floor(no/2);
            pf            = zeros(3,no);            
            pf(1,1:num)   = sqrt(2)/2*[0,linspace(0.25,0.5,num-1)];
            pf(1,num+1:no)= sqrt(2)/2*linspace(0.75,1.0,no-num);  
            pf(2,:)       = pf(1,:);   
            pf(3,:)       = 3*(1-(sqrt(2)*pf(1,:)).^2);  
        case 'MaOP4'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),[linspace(0,1/4,floor(num/2)),linspace(0.75,1,num-floor(num/2))]);
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = ps(1,:).*ps(2,:);
            pf(2,:)     = (1-ps(1,:)).*ps(2,:);
            pf(3,:)     = 1-ps(2,:);   
        case 'MaOP5'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid([linspace(0,1/4,floor(3*num/4)),linspace(3/4,1,num-floor(3*num/4))],linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
            pf(2,:)     = 4*cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
            pf(3,:)     = sin(0.5*pi*ps(1,:));   
        case 'MaOP6'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),[0,linspace(1/4,0.5,floor(num/2)),linspace(0.75,1,num-floor(num/2)-1)]);
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = ps(1,:).*ps(2,:);
            pf(2,:)     = (1-ps(1,:)).*ps(2,:);
            pf(3,:)     = 1-ps(2,:);   
        case 'MaOP7'
            num         = floor(sqrt(1000));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            tpf          = zeros(3,no);
            tpf(1,:)     = ps(1,:);
            tpf(2,:)     = ps(2,:);
            tpf(3,:)     = 3-sum(ps(1:2,:).*(1+sin(3*pi*ps(1:2,:))));
            pf           = n_sort(tpf);
            clear s t;
        case 'MaOP8'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),[linspace(0,1/4,floor(num/2)),linspace(0.75,1,num-floor(num/2))]);
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
            pf(2,:)     = cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
            pf(3,:)     = sin(0.5*pi*ps(1,:));  
        case {'gf10','MaOP9'}
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            ff          = zeros(3,no);
            ff(1,:)     = ps(1,:).*ps(2,:);
        for j=1:no
           if ps(1,j)>0.2&&ps(1,j)<=0.8&&ps(2,j)>0.2&&ps(2,j)<=0.8
              ff(2,j)   = (1.0 - ps(2,j))*ps(1,j) ;
              ff(3,j)   = 2.0-ps(1,j);
           else
              ff(2,j)   = (1.0 - ps(2,j))*ps(1,j) ;
              ff(3,j)   = 1.0-ps(1,j);
          end
        end
           pf           = n_sort(ff);
        case 'MaOP10'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            ff          = zeros(3,no);           
        for j=1:no
           if ps(1,j)>0.2&&ps(1,j)<=0.8&&ps(2,j)>0.2&&ps(2,j)<=0.8
              ff(1,j)   = 2*(1-cos(0.5*pi*ps(1,j))).*(1-cos(0.5*pi*ps(2,j)));
              ff(2,j)   = 2*(1-cos(0.5*pi*ps(1,j))).*(1-sin(0.5*pi*ps(2,j)));
              ff(3,j)   = 2*(1-sin(0.5*pi*ps(1,j)));
           else
              ff(1,j)   = (1-cos(0.5*pi*ps(1,j))).*(1-cos(0.5*pi*ps(2,j)));
              ff(2,j)   = (1-cos(0.5*pi*ps(1,j))).*(1-sin(0.5*pi*ps(2,j)));
              ff(3,j)   = (1-sin(0.5*pi*ps(1,j)));
           end
        end
           pf           = n_sort(ff);
            case 'gf13'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            ff          = zeros(3,no);           
        for j=1:no
           if ps(1,j)>0.4&&ps(1,j)<=0.6&&ps(2,j)>0.4&&ps(2,j)<=0.6
              ff(1,j)   = 2*cos(0.5*pi*ps(1,j)).*cos(0.5*pi*ps(2,j));
              ff(2,j)   = 2*cos(0.5*pi*ps(1,j)).*sin(0.5*pi*ps(2,j));
              ff(3,j)   = 2*sin(0.5*pi*ps(1,j));
           else
              ff(1,j)   = cos(0.5*pi*ps(1,j)).*cos(0.5*pi*ps(2,j));
              ff(2,j)   = cos(0.5*pi*ps(1,j)).*sin(0.5*pi*ps(2,j));
              ff(3,j)   = sin(0.5*pi*ps(1,j));
           end
        end
           pf           = n_sort(ff);
        case 'gf14'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid([linspace(0,0.25,floor(num/2)),linspace(0.75,1,num-floor(num/2))],linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = ps(1,:).*ps(2,:);
            pf(2,:)     = (1-ps(1,:)).*ps(2,:);
            pf(3,:)     = 1-ps(2,:);   
       case 'MaOP11'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            y           = zeros(3,no);            
            h           = zeros(1,no);
        for i=1:no
           if abs(ps(1,i)*ps(2,i)- ps(1,i)*(1- ps(2,i)))>0.01
              h(i)   =  max(0,(1+0.4)*sin(2*pi*(ps(1,i))));           
           end
        end    
        y(1,:)       = (1+h).*ps(1,:).*ps(2,:);
        y(2,:)       = (1+h).*ps(1,:).*(1-ps(2,:));
        y(3,:)       = (1+h).*(1-ps(1,:));
         pf          = n_sort(y);
       case 'MaOP12'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            y           = zeros(3,no);            
            h           = zeros(1,no);
        for i=1:no
           if abs(ps(1,i)*ps(2,i)- ps(1,i)*(1- ps(2,i)))>0.01
              h(i)   =  max(0,-(1+0.4)*sin(2*pi*(ps(1,i))));           
           end
        end    
        y(1,:)       = (1+h).*ps(1,:).*ps(2,:);
        y(2,:)       = (1+h).*ps(1,:).*(1-ps(2,:));
        y(3,:)       = (1+h).*(1-ps(1,:));
         pf           = n_sort(y);
       case 'MaOP13'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            y           = zeros(3,no);            
            h           = zeros(1,no);
        for i=1:no
           if abs(ps(1,i)*ps(2,i)- ps(1,i)*(1- ps(2,i)))>0.01
              h(i)   =  max(0,-(1+0.4)*sin(2*pi*(ps(1,i))));           
           end
        end    
        y(1,:)       = (1+h).*cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
        y(2,:)       = (1+h).*cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
        y(3,:)       = (1+h).*sin(0.5*pi*ps(1,:));
         pf           = n_sort(y);
        case 'MaOP14'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            y           = zeros(3,no);            
            h           = zeros(1,no);
        for i=1:no
           if abs(cos(0.5*pi*ps(1,i))*cos(0.5*pi*ps(2,i))- cos(0.5*pi*ps(1,i))*sin(0.5*pi*ps(2,i)))>0.01
              h(i)   =  max(0,(1+0.4)*sin(2*pi*(ps(1,i))));           
           end
        end    
        y(1,:)       = (1+h).*cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
        y(2,:)       = (1+h).*cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
        y(3,:)       = (1+h).*sin(0.5*pi*ps(1,:));
        pf           = n_sort(y);  
      case 'MaOP15'
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            y           = zeros(3,no);            
            h           = zeros(1,no);
        for i=1:no
           if abs(ps(1,i)*ps(2,i)- ps(1,i)*(1- ps(2,i)))>0.01
              h(i)   =  max(0,(1+0.4)*sin(2*pi*(ps(1,i))));           
           end
        end    
        y(1,:)       = (1+h).*cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:));
        y(2,:)       = (1+h).*cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:));
        y(3,:)       = (1+h).*sin(0.5*pi*ps(1,:));
        pf           = n_sort(y); 
    end
end
