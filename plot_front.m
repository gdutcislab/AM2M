% gudt EA lib by leichen
function plot_front(p,val,nrun)
PROBLEMS= ['MaOP1 '; 'MaOP2 '; 'MaOP3 '; 'MaOP4 '; 'MaOP5 '; 'MaOP6 '; 'MaOP7 '; 'MaOP8 '; 'MaOP9 '; 'MaOP10';'MaOP11';'MaOP12';'MaOP13';'MaOP14';'MaOP15';];
DIMX    = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
NOP     = [100 100 100 500 500 500 500 500 500 500 500 500 500 500 500 500]; 
PROPERTY= ['r.'; 'r.'; 'r.'; 'r.'; 'r.'; 'r.';'r.';'r.';'r.';'r.';'r.'; 'r.'; 'r.'; 'r.'; 'r.'; 'r.';];

%% 参数设置
    PF      = pareto( deblank(PROBLEMS(p,:)), NOP(p), DIMX(p) );    
    set(0,'units','centimeters');
    position=[10 10 6.5 6.5];
    h=figure;
    set(h,'PaperType','A4'); 
    set(h,'PaperUnits','centimeters'); 
    set(h,'paperpositionmode','auto');
    set(h,'PaperPosition',position);
    set(h,'units','centimeters');
    set(h,'position',position);
    hold off;
%% 删除劣解
       [dim,num_p]  = size(val);
        rank             = zeros(1,num_p);
        for k=1:num_p
            for j=k+1:num_p
                    if all(val(:,j)<=val(:,k))
                        rank(k)  = rank(k)+1;
                    elseif all(val(:,k)<=val(:,j))
                        rank(j)  = rank(j)+1;
                    end
            end
        end 
        val(:,rank~=0)=[];
%%画图
    if size(PF,1)== 3
        plot3(PF(1,:),PF(2,:),PF(3,:),deblank(PROPERTY(p,:)),'MarkerSize',1); hold on;
        scatter3(val(1,:),val(2,:),val(3,:),4);
    else
        plot(PF(1,:),PF(2,:),deblank(PROPERTY(p,:)),'MarkerSize',1); hold on;
        scatter(val(1,:),val(2,:),4);
    end
        
    set(gca,'FontSize',8);
    xlabel('f1');ylabel('f2');
    title(PROBLEMS(p,:));      
%     set(gca,'YTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
%     set(gca,'YTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});             
    if p<4
        set(gca,'ZTick',     [0 0.5 1.0 1.5 2.0 2.5 3.0]);
        set(gca,'ZTickLabel',{'0.0','0.5','1.0','1,5','2.0','2.5','3.0'});  
        view([122 20]);
        axis([0 1.2 0 1.2 0 3.2])
    elseif p==5
        xlim([0 4]); ylim([0 1]);  
        set(gca,'XTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
        set(gca,'XTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'}); 
        set(gca,'YTick',     [0 0.6 1.2 1.8 2.4 3.0 3.6 4.2]);
        set(gca,'YTickLabel',{'0.0','0.6','1,2','1.8','2.4','3.0','3.6','4.2'});
        set(gca,'ZTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
        set(gca,'ZTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});  
        view([134 12]);
        axis([0 1.2 0 4.2 0 1.2])
    elseif p==7
        set(gca,'ZTick',     [-1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6]);
        set(gca,'ZTickLabel',{'-1.2', '-0.8', '-0.4','0.0','0.4','0.8','1.2','1.6','2.0','2.4','2.8','3.2','3.6'});
        view([134 12]);
        axis([0 1.2 0 1.2 -1.2 3.6])
    else
        set(gca,'ZTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
        set(gca,'ZTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});  
        view([134 12]);
        axis([0 1.2 0 1.2 0 1.2])
    end
    grid off; box on;

%     subplot(1,2,2);
%     plot3(PS(1,:),PS(2,:),PS(3,:),deblank(PROPERTY(p,:)),'MarkerSize',4); hold on;  
%     set(gca,'FontSize',8);
%     xlabel('x1');ylabel('x2'); zlabel('x3');
%     title('Pareto set');
%     xlim(RA(1,:)); ylim(RA(2,:)); zlim(RA(3,:));            
%     view([-45 10]); grid off; box on;
    name=strcat(deblank(PROBLEMS(p,:)),'_R',num2str(nrun)); 
    f_1 = sprintf('pf_figs/%s.fig',name);
       saveas(h,f_1);              
    f_2 = sprintf('pf_figs/%s.eps',name);
       saveas(h,f_2);              

    close(h);
    set(0,'units','pixel');        
end
