function [params,mop,pop,inds,state]=inputparams(name_f,params,state)
    
    % 1）初始化真是的有效界面
    if( strcmpi(params.isCauIGD, 'yes'))
        PFStar   = load(strcat('PFStar/pf_',name_f,'.dat'));
        params.PFStar   =  PFStar';
    end
    
    % 2)初始化种群规模，最大进化代数，子种群的数目
    mop           = testmop(name_f);
    [params,mop]  = inputparameter(params,mop);    
    
    % 4)初始化权重
    if( strcmpi(params.useWeight, 'yes'))
        val_w                  = reference_point(mop);
        val_w                  = val_w./repmat(sqrt(sum(val_w.^2)),[mop.od,1]);
        params.weight          = val_w;
        params.popsize         = size(val_w,2);
    end
    
    % 5)初始化种群
    num_nod               = params.popsize;
    inds                  = randompoint(mop, num_nod);

    
    %% 评价每一个个体
    [inds,state]          = evaluate(inds,mop,state);    
    % 6)初始化中心点和子种群的中心规模
    val_cp                = Weight(params.num_class,mop.od);
    center                = val_cp./repmat(sqrt(sum(val_cp.^2)),[mop.od,1]);
    params.num_class      = size(center,2);
    sub                   = get_structure('subclass');
    pop                   = repmat(sub, [1,params.num_class]);
    for i=1:params.num_class
        pop(i).center     = center(:,i);
    end

    if( strcmpi(params.useWeight, 'yes'))
        team              = group(params,pop,val_w);
        for i=1:params.num_class
            pop(i).weight = 1./val_w(:,team{i});
            pop(i).num_ind= length(team{i});
        end
    else
        temv1     = mod(params.popsize,params.num_class);
        temv2     = floor(params.popsize/params.num_class);
        for i=1:temv1
            pop(i).num_ind=temv2+1;
        end
        for i=temv1+1:params.num_class
            pop(i).num_ind=temv2;
        end
    end   
    params.pmuta           = 1/mop.pd;
    params.idealmin        = 1000000000000*ones(mop.od,1);    
end
function [params,mop]=inputparameter(params,mop)
    switch  mop.name       
        case {'MaOP1','MaOP2','MaOP3'}
            params.num_class       = 10;
            params.iteration       = 110000;
            params.popsize         = 110;
        case{'MaOP4','MaOP5','MaOP6','MaOP7','MaOP8','MaOP9','MaOP10','MaOP11','MaOP12','MaOP13','MaOP14','MaOP15'}
            params.num_class       = 20;
            params.iteration       = 275*2000;
            params.popsize         = 275;
        case {'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'}
            params.num_class       = 10;
            params.iteration       = 210*750;
            params.popsize         = 210;
        case{'IN1','IN2','IN3','IN4','IN5','IN6'}
            params.num_class       = 17;
            params.iteration       = 600000;
            params.popsize         = 300;            
    end 
end

%% 产生个体
function ind = randompoint(prob, n)
    if (nargin==1)
        n=1;
    end
    randarray = rand(prob.pd, n);
    lowend = prob.domain(:,1);
    span = prob.domain(:,2)-lowend;
    point = randarray.*(span(:,ones(1, n)))+ lowend(:,ones(1,n));
    cellpoints = num2cell(point, 1);
    indiv = get_structure('individual');
    ind = repmat(indiv, [1, n]);
    [ind.parameter] = cellpoints{:};
end
function val_w = Weight(popsize,objDim)
    if objDim==2
        start      = 1/(popsize*100000);
        val_w(1,:) = linspace(start,1-start,popsize);
        val_w(2,:) = ones(1,popsize)-val_w(1,:);
    else
        val_w = lhsdesign(popsize, objDim, 'criterion','maximin','iterations', 1000)';    
    end
end
