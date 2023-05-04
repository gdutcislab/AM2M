function weight = GW(val,num)
     [dim,siz]         = size(val);
     noval             = val./repmat(sqrt(sum(val.^2)),[dim,1]);
     dis               = noval'*noval;
     weight            = zeros(dim,num);
     loc               = floor(rand*siz)+1;
     weight(:,1)       = noval(:,loc);
     pdis              = dis(loc,:);
    for i=2:num
        [ser,loc]      = min(pdis);
        weight(:,i)    = noval(:,loc);
        pdis           = max([pdis;dis(loc,:)]);
    end
end