function str = get_structure( name )
switch name 
    % state: optimization state of one generation
    case 'state'
    str = struct(...
          'currentGen', 1,...         % current generation number
          'reg_obj', 0,...            % number of objective function evaluation
          'totalTime', 0,...          % total time from the beginning
          'stopCriterion',1,...       % 终止条件
          'val_IGD',[],...            % 计算的IGD的值
          'val_H',[],...              % 当前的H的值
          'archive',[]);              % 精英集              
      
    case 'individual' 
        str = struct(...
            ...%%%%%%%%%%%%%%%%%%%%%%%%%%%最简形式的参数，
            'parameter',[],...            %parameter是自变量，
            'objective',[],...            %objective是目标函数 
            ...%%%%%%%%%%%%%%%%%%%%%%%%%%%私有参数，针对特定算法的
            ...%%%%%%%%%%%%%%%%%基于支配关系的算法
            'rank',[],...                 %个体的rank，非支配排序中用到
            'distance',[],...             %拥挤距离，非支配排序算法中用到
            ...%%%%%%%%%%%%%%%%基于指标的算法
            'inValue',[]);                %指标值           

    case 'parameter'
        str = struct(...
            ...%%%%%%%%%%%%%%%%%%%%%%%%%%%公有参数
            'popsize',[],...              % popsize是种群规模
            'isCauIGD','no',...           % 是否计算IGD
            'isDebug',1,...               % 是否是调试版本，1是0否
            'resultOut','no',...          % 结果的输出方式，是保存还是显示
            'iteration',[],...            % iteration最大目标函数的计算次数
            'PFStar',[],...               % 有效界面上真实的均匀分布的点
            'isCauH','yes',...            % 是否计算H度量值 
            'useArchive','no',...         % 是否对精英集进行存档
            ...%%%%%%%%%%%%%%%%%%%%%%%%%%%私有参数
            ...%%%%%%%%%%%%分区域和并行计算的参数
            'num_class',[],...            % num_class子区域的数目
            ...%%%%%%%%%%%%杂交变异种的参数
            'pmuta',[],...                % pmuta变异的概率，
            'delta',[],...                % delta杂交变异时的模拟退火因子，
            'selectPro',0.7,...           % 怎么选择参与杂交变异的个体，在分区域算法中用到 
            ...%%%%%%%%%%%%选择算子用到的参数
            'selectMethod','MM',...       % 选择方式：‘NS’：nsga-II的选择方式；
            ...                           %'WS'加权和选择方式；‘MM’极大极小选择方式
            ...                           % 'VEGA’：vega的选择方式
            ...%%%%%%%%%%%%在基于权重选择算法中用到的权重参数
            'idealmin',[],...              % 到目前为止所发现的最小值
            'nadirpoint',[],...            % nadirpoint
            'useWeight','yes',...          % 是否需要权重
            'normalization','yes',...      % 是否规范化
            'weight',[],...                % 所对应的权重
            ...%%%%%%%%%%%在hype中用到的参数
            'reference',[]);               % 计算H-value时用到的参考点          
            
    case 'testmop'
         str = struct(...
             'name',[],...                %name是测试函数的名称，
             'od',[],...                  %od是目标函数的维数，
             'pd',[],...                  %pd是自变量的维数，
             'domain',[],...              %domain是自变量的取值范围
             'func',[]);                  %func是测试函数
                                        

    case 'subclass' 
        str = struct(...
            ...%%%%%%%%%%%%%%%%%%%%%%%%%%%公有参数
            'center',[],...               %center是子区域的中心
            'inter',[],...                %inter是子种群
            'num_ind',[],...              % num_ind子种群的规模
            ...%%%%%%%%%%%%%%%%%%%%%%%%%%%私有参数
            ...%%%%%%%%%%%%%%采用权重设计选择算子的权重参数
            'weight',[]);               % 每个子种群的权重
      
    otherwise
        error('the structure name requried does not exist!');

end