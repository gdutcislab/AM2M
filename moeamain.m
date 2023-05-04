%% 更新日志：2017-2-20 更新 2017年5月8号测试WFG
clear; close; clc;                            
% name_func  = {'WFG2'};
name_func  = { 'MaOP1','MaOP2','MaOP3','MaOP4','MaOP5','MaOP6','MaOP7','MaOP8','MaOP9','MaOP10','MaOP11','MaOP12','MaOP13','MaOP14','MaOP15'};
   for nrun=1:1
        for seq = 1:1      
            name_f   = char(name_func(seq));
            seed     = 60+nrun;randn('state',seed);rand('state',seed);
            state    = get_structure( 'state');
            params   = get_structure( 'parameter');
            [params,mop,pop,inds,state]=inputparams(name_f,params,state);
            [params,mop,pop]=initialize(params,mop,pop,inds); 
            while state.stopCriterion 
                [params,mop,pop,state]=evolution(params,mop,pop,state);         
                fprintf('gen =   %d\n',state.currentGen);           
                state.currentGen=state.currentGen+1; 
                state  = checkstop(params,state);
            end        
            state=stateOutput(state,params,pop,mop,nrun);
        end
    end




  