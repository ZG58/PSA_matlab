function InputParams = V3_ProcessInputParameters(process_vars, material, N)
%InputParameters: 指定PSA模拟的输入参数
%   有限体积数始终是调用脚本/
%   函数的输入。此外，x 用于允许从调用脚本/
%   函数改变特定参数（时间、压力、吸附热等）。
%   对于简单模拟而言，这不是必需的。
    
%% 声明模拟的所有输入参数
    
    % 检索作为函数输入提供的所需过程变量
    L           = process_vars(1)       ;   % 吸附塔长度 [m]   
    P_0         = process_vars(2)       ;   % 吸附压力 [Pa]
    ndot_0      = process_vars(3)       ;   % 入口摩尔通量 [mol/s/m^2]
    t_ads       = process_vars(4)       ;   % 吸附步骤时间 [s]
    alpha       = process_vars(5)       ;   % 轻组分产品回流比 [-]
    beta        = process_vars(6)       ;   % 重组分产品回流比 [-]
    P_I         = process_vars(7)       ;   % 中间压力 [Pa]
	P_l         = process_vars(8)       ;   % 吹扫压力 [Pa]
    t_pres      = process_vars(9)       ;   % 加压步骤时间 [s]
    t_CnCdepres = process_vars(10)      ;   % 逆流降压步骤时间 [s]
    t_CoCdepres = process_vars(11)      ;   % 并流降压步骤时间 [s]
    
    % 检索取决于吸附剂材料的参数
    material_propertry = material{1}    ;
    IsothermPar        = material{2}    ; 
    
    % 操作床层参数
    % t_pres      = 20                    ;   % 最大/加压步骤时间 [s]
    % t_CnCdepres = 30                    ;   % 最大/逆流降压步骤时间 [s]
    % t_CoCdepres = 70                    ;   % 最大/并流降压步骤时间 [s]
    t_LR        = t_ads                 ;   % 轻组分回流步骤时间 [s]
    t_HR        = t_LR                  ;   % 重组分回流步骤时间 [s]
    tau         = 0.5                   ;   % 用于确定压力变化速度的参数
    P_inlet     = 1.02                  ;   % 吸附步骤入口处的进料气压力
    
    % 烟气参数和常数
    R          = 8.314                  ;   % 通用气体常数 [J/mol/K : Pa*m^3/mol/K]
    T_0        = 313.15                 ;   % 烟气进料温度 [K]
    y_0        = 0.15                   ;   % 入口气体CO2摩尔分数[-]
    Ctot_0     = P_0/R/T_0              ;   % 入口总浓度 [mol/m^3]
    v_0        = ndot_0/Ctot_0          ;   % 入口速度和比例参数 [m/s]
    mu         = 1.72e-5                ;   % 气体粘度 [Pa*s]
    epsilon    = 0.37                   ;   % 空隙率
    D_m        = 1.2995e-5              ;   % 分子扩散系数 [m^2/s]
    K_z        = 0.09                   ;   % 气相热导率 [W/m/k]
    C_pg       = 30.7                   ;   % 气体比热 [J/mol/k]
    C_pa       = 30.7                   ;   % 吸附相比热 [J/mol/k]
    MW_CO2     = 0.04402                ;   % CO2分子量 [kg/mol]
    MW_N2      = 0.02802                ;   % N2分子量 [kg/mol]
	%feed_gas  = 'Constant Pressure'    ;   % 进料步骤中的烟气是恒压还是恒速
    feed_gas   = 'Constant Velocity'    ;   % 进料步骤中的烟气是恒压还是恒速
    
    % 吸附剂参数
    % ro_s        = 1130                  ;   % 吸附剂密度 [kg/m^3]
    ro_s        = material_propertry(1) ;
    r_p         = 1e-3                  ;   % 颗粒半径 [m]
    C_ps        = 1070                  ;   % 吸附剂比热容 [J/kg/K]
    q_s         = 5.84                  ;   % 摩尔负载量比例因子 [mol/kg]
    q_s0        = q_s*ro_s              ;   % 摩尔负载量比例因子 [mol/m^3]
    k_CO2_LDF   = 0.1631                ;   % CO2的传质系数 [1/s]
    k_N2_LDF    = 0.2044                ;   % N2的传质系数 [1/s]
    
    % 等温线参数
    q_s_b      = [IsothermPar(1),  IsothermPar(7)]   ;   % b位上的饱和负载量 [mol/kg]
    q_s_d      = [IsothermPar(2),  IsothermPar(8)]   ;   % d位上的饱和负载量 [mol/kg]
    b          = [IsothermPar(3),  IsothermPar(9)]   ;   % b位的指前因子 [Pa-1]
    d          = [IsothermPar(4),  IsothermPar(10)]  ;   % d位的指前因子 [Pa-1]
    deltaU_b   = [IsothermPar(5),  IsothermPar(11)]  ;   % b位的吸附热 [J/mol]
    deltaU_d   = [IsothermPar(6),  IsothermPar(12)]  ;   % d位的吸附热 [J/mol]
    
    % deltaU     = [-36000, -15800]     ; 
    deltaU     = [material_propertry(2), material_propertry(3)]     ; 
    
    
    
%% 将值分配给必要的变量
    Params     = zeros(39, 1) ;
    Params(1)  = N			  ;
    Params(2)  = deltaU(1)    ;
    Params(3)  = deltaU(2)    ;
    Params(4)  = ro_s		  ;
    Params(5)  = T_0		  ;
    Params(6)  = epsilon	  ;
    Params(7)  = r_p		  ;
    Params(8)  = mu			  ;
    Params(9)  = R			  ;
    Params(10) = v_0		  ;
    Params(11) = q_s0		  ;
    Params(12) = C_pg		  ;
    Params(13) = C_pa		  ;
    Params(14) = C_ps		  ;
    Params(15) = D_m		  ;
    Params(16) = K_z		  ;
    Params(17) = P_0		  ;
    Params(18) = L			  ;
    Params(19) = MW_CO2		  ;
    Params(20) = MW_N2		  ;
    Params(21) = k_CO2_LDF	  ;
    Params(22) = k_N2_LDF	  ;
    Params(23) = y_0		  ;
    Params(24) = tau		  ;
    Params(25) = P_l		  ;
    Params(26) = P_inlet	  ;
    Params(27) = 1			  ;   % 吸附出口y的位置 = 轻组分回流入品y: y_LP
                                  % y_LR = 1 - 轻组分回流步骤中入口CO2摩尔分数的无初始猜测值
    Params(28) = 1			  ;   % 吸附出口T的位置 = 轻组分回流入品T: T_LP
                                  % T_LR = 1 - 轻组分回流步骤中入口温度的无初始猜测值
    Params(29) = 1			  ;   % 吸附出口ndot的位置 = 轻组分回流入品ndot
                                  % ndot_LR = 1 - 轻组分回流步骤中入口ndot的无初始猜测值
    Params(30) = alpha    	  ;
    Params(31) = beta         ;
    Params(32) = P_I          ;
    Params(33) = y_0          ;   % 逆流降压出口y的位置 = 重组分回流入品y: y_HP
                                  % y_HR = y_0 - 重组分回流步骤中入口CO2摩尔分数的初始猜测值
    Params(34) = T_0          ;   % 逆流降压出口T的位置 = 重组分回流入品T: T_HP
                                  % T_HR = T_0 - 重组分回流步骤中入口温度的初始猜测值
    Params(35) = ndot_0*beta  ;   % 300/30   % 逆流降压出口ndot的位置 = 重组分回流入品ndot
                                  % ndot_HR = ndot_0*beta - 重组分回流步骤中入口ndot的初始猜测值
    Params(36) = 0.01    	  ;   % 吸附出口y的位置 = 并流加压入口y: y_LP
                                  % y_LR = 0.01 - 并流加压步骤中入口CO2摩尔分数的初始猜测值
    Params(37) = T_0    	  ;   % 吸附出口T的位置 = 并流加压入口T: T_LP
                                  % T_LR = T_0 - 并流加压步骤中入口温度的初始猜测值
    Params(38) = ndot_0  	  ;   % 吸附出口ndot的位置 = 并流加压入口ndot
                                  % 注意: 未使用，似乎非必需。 ndot_LR = ndot_0 - 
                                  % 并流加压步骤中入口ndot的初始猜测值
    
    if strcmpi(feed_gas, 'Constant Pressure') == 1
        Params(end) = 1 ;
    elseif strcmpi(feed_gas, 'Constant Velocity') == 1
        Params(end) = 0 ;
    else
        error('请指定进料步骤的入口速度或压力是否恒定')
    end
	
    Times          = [ t_pres; t_ads; t_CnCdepres; t_LR; t_CoCdepres; t_HR ] ;
 
    IsothermParams = [q_s_b, q_s_d, b, d, deltaU_b, deltaU_d, IsothermPar(13)] ;
%   

%% 经济参数
    
    desired_flow                          = 100          ;   % 每塔烟气期望流速 [mol/s]
    electricity_cost                      = 0.07         ;   % 电价 [$/kWh]
    hour_to_year_conversion               = 8000         ;   % 一年中的总小时数，剩余时间假定为维护停机时间 [hr/year]
    life_span_equipment                   = 20           ;   % 除吸附剂外的所有设备寿命 [年]
    life_span_adsorbent                   = 5            ;   % 吸附剂寿命 [年]
    CEPCI                                 = 536.4        ;   % 当年月份的CEPCI（2016年1月）。
    
    % 根据要模拟的循环更改此项
    cycle_time = t_pres + t_ads + t_HR + t_CoCdepres + t_CnCdepres + t_LR ;   % 1个循环所需的总时间 [s]
    
    EconomicParams    = zeros(6, 1)              ;
    EconomicParams(1) = desired_flow             ;
    EconomicParams(2) = electricity_cost         ;
    EconomicParams(3) = cycle_time               ;
    EconomicParams(4) = hour_to_year_conversion  ;
    EconomicParams(5) = life_span_equipment      ;
    EconomicParams(6) = life_span_adsorbent      ;
    EconomicParams(7) = CEPCI                    ;
%   
%% 将所有列表合并到一个可以轻松传递的元胞变量中
    InputParams{1} = Params          ;
    InputParams{2} = IsothermParams  ;
    InputParams{3} = Times           ;
    InputParams{4} = EconomicParams  ;
%     InputParams{5} = economic_class  ;
%   
end