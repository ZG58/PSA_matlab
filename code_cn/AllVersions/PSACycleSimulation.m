function [ objectives, constraints ] = PSACycleSimulation( x, material, type, N)

    % % 检索作为函数输入提供的所需过程变量
    % L           = process_vars(1)       ;   % 吸附塔长度 [m]   
    % P_0         = process_vars(2)       ;   % 吸附压力 [Pa]
    % ndot_0      = process_vars(3)       ;   % 入口摩尔通量 [mol/s/m^2]
    % t_ads       = process_vars(4)       ;   % 吸附步骤时间 [s]
    % alpha       = process_vars(5)       ;   % 轻组分产品回流比 [-]
    % beta        = process_vars(6)       ;   % 重组分产品回流比 [-]
    % P_I         = process_vars(7)       ;   % 中间压力 [Pa]
	% P_l         = process_vars(8)       ;   % 吹扫压力 [Pa]

    % x(1)                                    % 吸附压力 [Pa]
    % x(2)                                    % 吸附步骤的时间 [s] 
    % x(3)                                    % 轻组分产品回流比 [-] 
    % x(4)                                    % 进料速度 [m/s] 
    % x(5)                                    % 重组分产品回流比 [-] 
    % x(6)                                    % 吹扫压力 [Pa]

    process_variables = [1.0, x(1), x(1)*x(4)/8.314/313.15,  x(2),  x(3),  x(5), 1e4, x(6)] ;
	
	try
    [objectives, constraints] = PSACycle(process_variables, material, [], type, N) ;
	catch
    %warning('函数使用出现问题。为目标和约束违反赋零值');
    objectives(1)  = 1e5 ;
	objectives(2)  = 1e5 ;
	constraints(1) = 1 ;
	constraints(2) = 1 ;
	constraints(3) = 1 ;
	end
	
end