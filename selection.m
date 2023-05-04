function[nextpop,pop]=selection(params,pop,oldpop)
     if(nargin < 3)
         num = params.popsize;   
     else
         num = oldpop.num_ind;
     end
     if ( strcmpi(params.selectMethod, 'NS'))
         pop            = ndsort(pop);
         pop            = selectNS(pop,num);     
     elseif ( strcmpi(params.selectMethod, 'WS'))
         % Todo      
     elseif ( strcmpi(params.selectMethod, 'MM'))
         [nextpop,pop] = selectMM(pop,params,oldpop);
     elseif ( strcmpi(params.selectMethod, 'VEGA'))
         % Todo 
     else
         error( 'No defined the select method');
     end
end
%% nondominated sorting base selection
function nextpop = selectNS(pop,num)
    nextpop = pop(1:num);   %just for initializing
    rankVector = [pop.rank];
    n = 0;                  % individuals number of next population
    rank = 1;               % current rank number
    idx = find(rankVector == rank);
    numInd = length(idx);   % number of individuals in current front
    while( n + numInd <= num )
        nextpop( n+1 : n+numInd ) = pop( idx );
        n = n + numInd;
        rank = rank + 1;
        idx = find(rankVector == rank);
        numInd = length(idx);
    end
    if( n < num )
        distance   = [pop(idx).distance];
        distance   = [distance', idx'];
        distance   = sortrows( distance, -1);      
        idxSelect  = distance( 1:num-n, 2);          
        nextpop(n+1 : num) = pop(idxSelect);
    end
end
%% Tchebycheff decomposition based selection
function [nextpop,pop]= selectMM(pop,params,oldpop)   
         val            = [pop.objective];
         num_nod        = size(val,2);
         nextpop        = pop(1:oldpop.num_ind);
     if (strcmpi(params.normalization, 'yes'))    
         temp           = params.nadirpoint-params.idealmin;
         pos            = (temp < 1.0E-5);
         temp(pos)      = 1.0E-5;
         noval          = (val-repmat(params.idealmin,[1,num_nod]))./(repmat(temp,[1,num_nod]));
     else
         noval          = val-repmat(params.idealmin,[1,num_nod]);
     end
     weight              = oldpop.weight;
    for j=1:oldpop.num_ind
        max_evo          = max(repmat(weight(:,j),[1,num_nod]).*noval);
        [cfval,cfloc]    = min(max_evo);
        nextpop(j)       = pop(cfloc);
        pop(cfloc)       = [];
        noval(:,cfloc)   = [];
        num_nod          = num_nod-1;
    end
end





