function [params,mop,pop]= initialize(params,mop,pop,inds)         
         v                 = [inds.objective];
         params.idealmin   = min([params.idealmin,v],[],2);
         params.nadirpoint = max(v,[],2); 
    %% population initialization
         if (strcmpi(params.normalization, 'yes'))         
              noval  = (v-repmat(params.idealmin,[1,params.popsize]))./(repmat(params.nadirpoint-params.idealmin,[1,params.popsize]));
         else
              noval  = v-repmat(params.idealmin,[1,params.popsize]);
         end
             team    = group(params,pop,noval);
         for i=1:params.num_class
             num_p             = length(team{i});
             if num_p<=pop(i).num_ind
                 tst                      = floor(params.popsize*rand(1,pop(i).num_ind-num_p))+1;
                 selind                   = inds([team{i},tst]);
                 pop(i).inter             = selind;
             else
                 selind                   = inds(team{i});
                 pop(i).inter             = selection(params,selind,pop(i));
             end
         end
    %% initialize the Archive
    if ( strcmpi(params.useArchive, 'yes'))
        nspop           = ndsort(params,[pop.inter]);
        state.archive   = nspop([nspop.rank]==1); 
    end 
end




