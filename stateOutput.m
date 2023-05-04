function state=stateOutput(state,params,pop,mop,nrun)
    if params.isDebug==1||~state.stopCriterion
        gen  = state.currentGen;
        if (gen>0&&gen<11)||(mod(gen,10)==0&&gen<100)||(mod(gen,100)==0)||~state.stopCriterion
        % 当前种群的自变量与目标函数值
            if isempty(state.archive)
                individual    = [pop.inter];
                valf          = [individual.objective];
                valx          = [individual.parameter];
            else
                individual    = [state.archive];
                valf          = [individual.objective];
                valx          = [individual.parameter];
            end  
            %输出当前状态
            output(valf,valx,nrun,params,mop,state);
            plot_front(seq,valf,nrun);
        end 
    end
end
function output(val_f,val_x,nrun,params,mop,state)
        name      = mop.name;  
%% 保存当前的和最后的种群自变量的值
    if ~state.stopCriterion ||  ( strcmpi(params.resultOut, 'save'))
        file      = val_x';
        filename  = strcat('PS_',name,'_R',num2str(nrun),'gen',num2str(state.currentGen)); 
        filetype  = 'dat';
        folder    = 'PS_data';
        mySave(file,filename,filetype,folder)
    end
%% 保存当前的和最后的种群自变量的值
    if ~state.stopCriterion ||  ( strcmpi(params.resultOut, 'save'))
        file      = val_f';
        filename  = strcat('PF_',name,'_R',num2str(nrun),'gen',num2str(state.currentGen)); 
        filetype  = 'dat';
        folder    = 'PF_data';
        mySave(file,filename,filetype,folder)
    end
end

function mySave(file,filename,filetype,folder)
    name      = strcat(filename,'.',filetype);
    if ( strcmpi(filetype, 'dat'))
        fid       = fopen(name,'w');
        [numline,numi]   = size(file);
        for i=1:numline
            for j=1:numi
                fprintf(fid,'%8.6f ',file(i,j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
    elseif  ( strcmpi(filetype, 'fig'))
        saveas(file,name);
    end
    judge     = exist(folder);
    if judge ~= 7
        system(['mkdir ', folder]);
    end
    file_path = strcat(cd,'\',folder);
    movefile(name,file_path); 
end