    %% 1. 模拟并流加压步骤
        [t1, a] = ode15s(CoCPressurization_fxn, [0 tau_CoCPres], x0, opts1) ;
        
        % 修正输出（清理模拟结果）
        idx             = find(a(:, 1) < a(:, 2))         ;  % P_1  < P_2
        a(idx ,1)       = a(idx, 2)                       ;  % P_1  = P_2
        a(idx, N+3)     = a(idx, N+4)                     ;  % y_1  = y_2
        a(idx, 4*N+9)   = a(idx, 4*N+10)                  ;  % T_1  = T_2
        a(:, 2*N+5)     = a(:, 2*N+6)                     ;  % x1_1 = x1_2
        a(:, 3*N+7)     = a(:, 3*N+8)                     ;  % x2_1 = x2_2
        a(:, 3*N+6)     = a(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        a(:, 4*N+8)     = a(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        a(:, N+3:2*N+4) = max(min(a(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        % 存储并流加压步骤的最终条件 - 所有迭代
        % 以及塔前、后端处的CO2和总摩尔数
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t1*L/v_0, a, 'HPEnd') ;
        [totalEnd, CO2End, ~]     = StreamCompositionCalculator(t1*L/v_0, a, 'LPEnd') ;
        a_fin = [a_fin; a(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;

        fprintf('加压步骤totalFront = %.6f, totalEnd = %.6f\n', totalFront, totalEnd);

        
        % 为吸附步骤准备初始条件
        x10         = a(end, :)' ;  % 上一步的最终状态是
                                    % 当前步骤的初始状态
        x10(1)      = P_inlet    ;  % BC z=0 P: P_1   = P_inlet
        x10(N+2)    = 1          ;  % BC z=1 P: P_N+2 = 1
        x10(N+3)    = y_0        ;  % BC z=0 y: y_1   = y_0
        x10(2*N+4)  = x10(2*N+3) ;  % BC z=1 y: y_N+2 = y_N+1
        x10(4*N+9)  = 1          ;  % BC z=0 T: T_1   = 1
        x10(5*N+10) = x10(5*N+9) ;  % BC z=1 T: T_N+2 = T_N+1
        
        % 存储吸附步骤的初始条件 - 所有迭代
        b_in = [b_in; x10'] ;
        
        % PSA循环第一步的状态初始条件
        statesIC = a(1, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9]) ;
%       
    %% 2. 模拟吸附步骤
        [t2, b] = ode15s(Adsorption_fxn, [0 tau_ads], x10, opts2) ;
        
        % 修正输出（清理模拟结果）
        idx             = find(b(:, N+1) < 1)             ;  % P_N+1 < 1 = P_N+2
        b(idx, N+2)     = b(idx, N+1)                     ;  % P_N+2 = P_N+1
        b(:, 2*N+5)     = b(:, 2*N+6)                     ;  % x1_1 = x1_2
        b(:, 3*N+7)     = b(:, 3*N+8)                     ;  % x2_1 = x2_2
        b(:, 3*N+6)     = b(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        b(:, 4*N+8)     = b(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        b(:, N+3:2*N+4) = max(min(b(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        if Params(end) == 0
            %b = VelocityCorrection(b, ndot_0, 'HPEnd') ;
            b = velocitycleanup(b)                     ;
        end
        
        % 存储吸附步骤的最终条件 - 所有迭代
        % 以及塔前、后端处的CO2和总摩尔数
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t2*L/v_0, b, 'HPEnd') ;
        [totalEnd, CO2End, TEnd]  = StreamCompositionCalculator(t2*L/v_0, b, 'LPEnd') ;
        b_fin = [b_fin; b(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;

        fprintf('吸附步骤totalFront = %.6f, totalEnd = %.6f\n', totalFront, totalEnd);

        
        % 添加并更新轻组分回流步骤所需的参数。这些
        % 是吸附步骤出口（轻组分产品端）的组成和温度，
        % 它们是轻组分回流步骤的入口参数
        y_LR        = CO2End/totalEnd ;
        T_LR        = TEnd            ;
        ndot_LR     = totalEnd/t_ads  ;
        Params(27)  = y_LR            ;
        Params(28)  = T_LR            ;
        Params(29)  = ndot_LR         ;
        
        % 使用更新后的参数调用轻组分回流步骤的函数
        LightReflux_fxn = @(t, x) FuncLightReflux(t, x, Params, IsothermParams) ;
        
        % 为重组分回流步骤准备初始条件
        x20         = b(end, :)' ;  % 上一步的最终状态是
                                    % 当前步骤的初始状态
        x20(1)      = P_inlet    ;  % BC z=0 P: P_1   = P_inlet
        x20(1)      = x20(2)     ;   
        x20(N+2)    = 1          ;  % BC z=1 P: P_N+2 = 1
        x20(N+3)    = y_HR1       ;  % BC z=0 y: y_1   = y_HR1
        x20(2*N+4)  = x20(2*N+3) ;  % BC z=1 y: y_N+2 = y_N+1
        x20(4*N+9)  = T_HR1/T_0   ;  % BC z=0 T: T_1   = T_HR1/T_0
        x20(5*N+10) = x20(5*N+9) ;  % BC z=1 T: T_N+2 = T_N+1
        
        % 存储重组分回流步骤的初始条件 - 所有迭代
        c_in = [c_in; x20'] ;
%       
    %% 3. 模拟重组分回流步骤
        [t3, c] = ode15s(HeavyReflux1_fxn, [0 tau_HR1], x20, opts3) ;
        
        % 修正输出（清理模拟结果）
        idx             = find(c(:, N+1) < 1)             ;  % P_N+1 < 1 = P_N+2
        c(idx, N+2)     = c(idx, N+1)                     ;  % P_N+2 = P_N+1
        c(:, 2*N+5)     = c(:, 2*N+6)                     ;  % x1_1 = x1_2
        c(:, 3*N+7)     = c(:, 3*N+8)                     ;  % x2_1 = x2_2
        c(:, 3*N+6)     = c(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        c(:, 4*N+8)     = c(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        c(:, N+3:2*N+4) = max(min(c(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        if Params(end) == 0
            c = VelocityCorrection(c, ndot_HR1, 'HPEnd') ;
            %c = velocitycleanup(c)                      ;
        end
        
        % 存储重组分回流步骤的最终条件 - 所有迭代
        % 以及塔前、后端处的CO2和总摩尔数
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t3*L/v_0, c, 'HPEnd') ;
        [totalEnd, CO2End, ~]     = StreamCompositionCalculator(t3*L/v_0, c, 'LPEnd') ;
        c_fin = [c_fin; c(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;

        fprintf('重回流 totalFront = %.6f, totalEnd = %.6f\n', totalFront, totalEnd);

        
        % 为并流降压步骤准备初始条件
        x30         = c(end,:)'   ;   % 上一步的最终状态是
                                      % 当前步骤的初始状态
        x30(1)      = x30(2)      ;   % BC z=0 P: P_1   = P_2
        x30(N+2)    = x30(N+1)    ;   % BC z=1 P: P_N+2 = P_N+1
        x30(N+3)    = x30(N+4)    ;   % BC z=0 y: y_1   = y_2
        x30(2*N+4)  = x30(2*N+3)  ;   % BC z=1 y: y_N+2 = y_N+1
        x30(4*N+9)  = x30(4*N+10) ;   % BC z=0 T: T_1   = T_2
        x30(5*N+10) = x30(5*N+9)  ;   % BC z=1 T: T_N+2 = T_N+1
        
        % 存储并流降压步骤的初始条件 - 所有迭代
        d_in = [d_in; x30'] ;
%
    %% 4. 模拟并流降压步骤
        [t4, d] = ode15s(CoCDepressurization_fxn, [0 tau_CoCDepres], x30, opts4) ;
        
        % 修正输出（清理模拟结果）
        idx             = find(d(:, N+1) < d(:, N+2))         ;  % P_N+1 < 1 = P_N+2
        d(idx ,N+2)     = d(idx, N+1)                     ;  % P_N+2 = P_N+1
        d(:, 2*N+5)     = d(:, 2*N+6)                     ;  % x1_1 = x1_2
        d(:, 3*N+7)     = d(:, 3*N+8)                     ;  % x2_1 = x2_2
        d(:, 3*N+6)     = d(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        d(:, 4*N+8)     = d(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        d(:, N+3:2*N+4) = max(min(d(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        % 存储逆流降压步骤的最终条件 - 所有迭代
        % 以及塔前、后端处的CO2和总摩尔数
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t4*L/v_0, d, 'HPEnd') ;
        [totalEnd, CO2End, ~]     = StreamCompositionCalculator(t4*L/v_0, d, 'LPEnd') ;
        d_fin = [d_fin; d(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;

        fprintf('并流降压totalFront = %.6f, totalEnd = %.6f\n', totalFront, totalEnd);

           
        % 为逆流降压步骤准备初始条件
        x50         = d(end,:)'   ;   % 上一步的最终状态是
                                      % 当前步骤的初始状态
        x50(1)      = x50(2)      ;   % BC z=0 P: P_1   = P_2
        x50(N+2)    = x50(N+1)    ;   % BC z=1 P: P_N+2 = P_N+1
        x50(N+3)    = x50(N+4)    ;   % BC z=0 y: y_1   = y_2
        x50(2*N+4)  = x50(2*N+3)  ;   % BC z=1 y: y_N+2 = y_N+1
        x50(4*N+9)  = x50(4*N+10) ;   % BC z=0 T: T_1   = T_2
        x50(5*N+10) = x50(5*N+9)  ;   % BC z=1 T: T_N+2 = T_N+1
        
        % 存储轻组分回流步骤的初始条件 - 所有迭代
        f_in = [f_in; x50'] ;
% 
    %% 5. 模拟逆流降压步骤
        [t6, f] = ode15s(CnCDepressurization_fxn, [0 tau_CnCDepres], x50, opts6) ;
        
        % 修正输出（清理模拟结果）
        idx             = find(f(:, 2) < f(:, 1))         ;  % P_2  < P_1
        f(idx ,1)       = f(idx, 2)                       ;  % P_1  = P_2
        f(:, 2*N+5)     = f(:, 2*N+6)                     ;  % x1_1 = x1_2
        f(:, 3*N+7)     = f(:, 3*N+8)                     ;  % x2_1 = x2_2
        f(:, 3*N+6)     = f(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        f(:, 4*N+8)     = f(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        f(:, N+3:2*N+4) = max(min(f(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        % 存储逆流降压步骤的最终条件 - 所有迭代
        % 以及塔前、后端处的CO2和总摩尔数
        [totalFront, CO2Front, TFront] = StreamCompositionCalculator(t6*L/v_0, f, 'HPEnd') ;
        [totalEnd, CO2End, ~]     = StreamCompositionCalculator(t6*L/v_0, f, 'LPEnd') ;
        f_fin = [f_fin; f(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;

        fprintf('逆流降压totalFront = %.6f, totalEnd = %.6f\n', totalFront, totalEnd);


        % 计算重组分回流步骤所需的参数
        y_HR1       = CO2Front/totalFront    ;
        T_HR1       = TFront                 ;
        ndot_HR1    = totalFront.*beta/t_HR1  ;
        Params(33) = y_HR1                   ;
        Params(34) = T_HR1                   ;
        Params(35) = ndot_HR1                ;
        
        HeavyReflux1_fxn = @(t, x) FuncHeavyReflux1(t, x, Params, IsothermParams) ;
        
        % 为轻组分回流步骤准备初始条件
        x60         = f(end,:)'    ;  % 上一步的最终状态是
                                      % 当前步骤的初始状态
        x60(1)      = P_l/P_0      ;  % BC z=0 P: P_1   = P_l/P_0
        x60(N+3)    = x60(N+4)     ;  % BC z=0 y: y_1   = y_2
        x60(2*N+4)  = y_LR         ;  % BC z=1 y: y_N+2 = y_LR
        x60(4*N+9)  = x60(4*N+10)  ;  % BC z=0 T: T_1   = T_2
        x60(5*N+10) = T_LR/T_0     ;  % BC z=1 T: T_N+2 = T_LR/T_0
        
        % 存储轻组分回流步骤的初始条件 - 所有迭代
        g_in = [g_in; x60'] ;
%       
    %% 6. 模拟轻组分回流步骤
        [t7, g] = ode15s(LightReflux_fxn, [0 tau_LR], x60, opts7) ;
        
        % 修正输出（清理模拟结果）
        idx             = find(g(:, 2) < g(:, 1))         ;  % P_2  < P_1
        g(idx ,1)       = g(idx, 2)                       ;  % P_1  = P_2
        g(:, 2*N+5)     = g(:, 2*N+6)                     ;  % x1_1 = x1_2
        g(:, 3*N+7)     = g(:, 3*N+8)                     ;  % x2_1 = x2_2
        g(:, 3*N+6)     = g(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        g(:, 4*N+8)     = g(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        g(:, N+3:2*N+4) = max(min(g(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        g = VelocityCorrection(g, ndot_LR*alpha, 'LPEnd') ;                        
        
        % 存储轻组分回流步骤的最终条件 - 所有迭代
        % 以及塔前、后端处的CO2和总摩尔数
        [totalFront, CO2Front, TFront]  = StreamCompositionCalculator(t7*L/v_0, g, 'HPEnd') ;
        [totalEnd, CO2End, ~]           = StreamCompositionCalculator(t7*L/v_0, g, 'LPEnd') ;
        g_fin = [g_fin; g(end, :), CO2Front, totalFront, CO2End, totalEnd]                  ;

        fprintf('轻回流 totalFront = %.6f, totalEnd = %.6f\n', totalFront, totalEnd);

        
        % 为并流加压步骤准备初始条件
        x0         = g(end, :)' ;  % 上一步的最终状态是
                                   % 当前步骤的初始状态
        x0(1)      = x0(2)      ;  % BC z=0 P: P_1   = P_2
        x0(N+2)    = x0(N+1)    ;  % BC z=1 P: P_N+2 = P_N+1
        x0(N+3)    = y_0        ;  % BC z=0 y: y_1   = y_0
        x0(2*N+4)  = x0(2*N+3)  ;  % BC z=1 y: y_N+2 = y_N+1
        x0(4*N+9)  = 1          ;  % BC z=0 T: T_1   = 1
        x0(5*N+10) = x0(5*N+9)  ;  % BC z=1 T: T_N+2 = T_N+1
        
        % PSA循环最后一步的状态最终条件
        statesFC = g(end, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9]) ;
%