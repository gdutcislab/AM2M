function [params,mop,pop,state]= evolution(params,mop,pop,state)

    %% 1)crossover & mutation
         params.delta          = (1- state.reg_obj/params.iteration)^0.7;
         [newpop,state]        = CMOp(params,mop,pop,state);
    
    %% 2)update population by selection
        if mod(state.reg_obj/params.popsize,100)==0
           [pop,params,state]  = dextractPop(pop,newpop,params,state);
        else
           [pop,params,state]  = extractPop(pop,newpop,params,state);
       end    
end
%% general update
function [pop,params,state] = extractPop(pop,newpop,params,state) 
     allind              = [pop.inter,newpop.inter];
     val                 = [allind.objective];
     num_nod             = size(val,2);
     params.idealmin     = min([params.idealmin,val],[],2);
     params.nadirpoint   = max(val,[],2); 
     if (strcmpi(params.normalization, 'yes')) 
          temp           = params.nadirpoint-params.idealmin;
          pos            = (temp < 1.0E-5);
          temp(pos)      = 1.0E-5;
          noval          = (val-repmat(params.idealmin,[1,num_nod]))./(repmat(temp,[1,num_nod]));
     else
         noval           = val-repmat(params.idealmin,[1,num_nod]);
     end
     team                = group(params, pop, noval);
     for i=1:params.num_class         
            num_p             = length(team{i});
           if num_p<=pop(i).num_ind
               tst            = floor(num_nod*rand(1,pop(i).num_ind-num_p))+1;
               selind         = allind([team{i},tst]);
               pop(i).inter   = selind;
           else
               selind         = allind(team{i});
               pop(i).inter   = selection(params,selind,pop(i));
           end         
     end
end
%% adaptive update
function [pop,params,state] = dextractPop(pop,newpop,params,state)   
         allind             = [pop.inter,newpop.inter];
         val                = [allind.objective];
         num_nod            = size(val,2);
         params.idealmin    = min([params.idealmin,val],[],2);
         params.nadirpoint  = max(val,[],2);
         if (strcmpi(params.normalization, 'yes'))
             temp           = params.nadirpoint-params.idealmin;
             pos            = (temp < 1.0E-5);
             temp(pos)      = 1.0E-5;
             noval          = (val-repmat(params.idealmin,[1,num_nod]))./(repmat(temp,[1,num_nod]));
         else
             noval          = val-repmat(params.idealmin,[1,num_nod]);
         end
        %% adaptive direction vector design
        weight                    = GW(noval,params.popsize);
        for i=1:params.num_class
           pop(i).center          = weight(:,i+1);
        end
        team                      = group(params,pop,noval);
        tm                        = zeros(1,params.num_class);
       %% identify the number of weights in each subregion
        for i=1:params.num_class
             tm(i)=length(team{i});
             pop(i).num_ind=floor(length(team{i})/2); 
         end
          if sum([pop.num_ind])<params.popsize
             temp        = params.popsize-sum([pop.num_ind]);
             [res,index] = sort(tm);
             t_index     = 0;           
                for i=1:params.num_class
                   if mod(tm(index(i)),2)==1&&t_index<temp
                       pop(index(i)).num_ind=pop(index(i)).num_ind+1;
                       t_index=t_index+1;
                   end
                end           
          end     
       %% adaptive weights setting&selection            
       for i=1:params.num_class  
           weight                = GW(noval(:,team{i}),pop(i).num_ind);
           pos                   = (weight < 1.0E-5);
           weight(pos)           = 1.0E-5;
           pop(i).weight         = 1./weight;           
           pop(i).inter          = selection(params,allind,pop(i)); 
       end        
end

%% generation of offsprings
function [newpop,state] = CMOp(params,mop,pop,state)
    PM          = [pop.inter];
    newpop      = pop;
    x_min    = mop.domain(:,1);
    x_max    = mop.domain(:,2);
    for i=1:params.num_class
        newinter  = pop(i).inter;
        for j=1:pop(i).num_ind
            parent1=pop(i).inter(j);
            if rand>params.selectPro
                parent2  = PM(floor(rand*params.popsize)+1);
            else
                loc                 = floor(rand*pop(i).num_ind)+1;
                parent2             = pop(i).inter(loc);
            end
            newinter(j).parameter   = SBXCM(parent1.parameter,parent2.parameter,x_min,x_max,1/mop.pd);
        end
        [newinter,state]            = evaluate(newinter,mop,state);
        newpop(i).inter             = newinter;
    end
end
%%