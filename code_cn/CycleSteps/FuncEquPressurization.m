function derivatives = FuncEquPressurization(~, state_vars, params, isotherm_params)
%#codegen
%均压升：计算均压升步骤中状态变量的变化
%   与MATLAB ode求解器（ode15s）结合使用，求解变压吸附（PSA）步骤的偏微分方程组。
%   此步骤中，塔底关闭，气体从塔顶进入。
%   该偏微分方程组通过使用有限体积法（FVM）求解，
%   其中空间域被离散为N个体积单元。
%   假设每个体积单元内的状态变量值是均匀的。
%   加权基本无振荡（WENO）格式用于通过计算有限体积单元壁面上的
%   状态变量值来提高计算精度。
%   
%   输入:
%       ~: 时间变量。由于在计算导数时未使用时间，
%          但ode求解器需要一个占位符，因此使用此符号。
%   
%       state_vars: 正在计算的当前时间的无量纲状态变量。
%          顺序如下 [P_1, P_2,...,P_N+2, y_1,
%          ...,y_N+2, q_CO2_1,...,q_CO2_N+2, q_N2_1,...,q_N2_N+2, T_1,...,
%          T_N+2]
%   
%       params: 模拟所需的所有参数。在
%          ProcessInputParameters函数文件中提供。
%   
%       isotherm_params: 等温线所需的所有参数。在
%          ProcessInputParameters函数文件中提供。
%   
%   输出:
%       derivatives: 状态变量的时间导数。
%       其顺序与状态变量相同。
%   
%% 检索过程参数
    N				=	params(1)	;
    deltaU_1	    =	params(2)	;
    deltaU_2	    =	params(3)	;
    ro_s			=	params(4)	;
    T_0				=	params(5)	;
    epsilon			=	params(6)	;
    r_p				=	params(7)	;
    mu				=	params(8)	;
    R				=	params(9)	;
    v_0				=	params(10)	;
    q_s0			=	params(11)	;
    C_pg			=	params(12)	;
    C_pa			=	params(13)	;
    C_ps			=	params(14)	;
    D_m				=	params(15)	;
    K_z				=	params(16)	;
    P_0				=	params(17)	;
    L				=	params(18)	;
    MW_CO2			=	params(19)	;
    MW_N2			=	params(20)	;
    k_1_LDF		    =	params(21)	;
    k_2_LDF		    =	params(22)	;
    tau				=	params(24)	;
    P_eq			=	params(36)	; % 均压目标压力
    y_eq            =   params(37)  ; % 均压进气组分
    T_eq            =   params(38)  ; % 均压进气温度
%   
%% 初始化状态变量
    P  = zeros(N+2, 1) ;
    y  = zeros(N+2, 1) ;
    x1 = zeros(N+2, 1) ;
    x2 = zeros(N+2, 1) ;
    T  = zeros(N+2, 1) ;
    
    P(1:N+2)  = state_vars(1:N+2)               ;
    y(1:N+2)  = max(state_vars(N+3:2*N+4), 0)   ;
    x1(1:N+2) = max(state_vars(2*N+5:3*N+6), 0) ;
    x2(1:N+2) = state_vars(3*N+7:4*N+8)         ;
    T(1:N+2)  = state_vars(4*N+9:5*N+10)        ;
%   
%% 初始化函数中使用的所有变量
    % 时间导数
    derivatives = zeros(5*N+10, 1) ;
    dPdt        = zeros(N+2, 1)    ;
    dPdt1       = zeros(N+2, 1)    ;
    dPdt2       = zeros(N+2, 1)    ;
    dPdt3       = zeros(N+2, 1)    ;
    dydt        = zeros(N+2, 1)    ;
    dydt1       = zeros(N+2, 1)    ;
    dydt2       = zeros(N+2, 1)    ;
    dydt3       = zeros(N+2, 1)    ;
    dx1dt       = zeros(N+2, 1)    ;
    dx2dt       = zeros(N+2, 1)    ;
    dTdt        = zeros(N+2, 1)    ;
    dTdt1       = zeros(N+2, 1)    ;
    dTdt2       = zeros(N+2, 1)    ;
    dTdt3       = zeros(N+2, 1)    ;
    % 空间导数
    dpdz        = zeros(N+2, 1)    ;
    dpdzh       = zeros(N+1, 1)    ;
    dydz        = zeros(N+2, 1)    ;
    d2ydz2      = zeros(N+2, 1)    ;
    dTdz        = zeros(N+2, 1)    ;
    d2Tdz2      = zeros(N+2, 1)    ;
%   
%% 计算所有使用的参数
    dz   = 1/N                                ;
    D_l  = 0.7*D_m + v_0*r_p                  ;
    Pe   = v_0*L/D_l                          ;
    phi  = R*T_0*q_s0*(1-epsilon)/epsilon/P_0 ;
    ro_g = P(1:N+2).*P_0/R./T(1:N+2)/T_0      ;
%   
%% 边界条件
%   均压升的边界条件，对应于塔底
%   （重组分产品端）关闭而塔顶（轻组分产品端）
%   打开以补充气体的情况。
%   对于塔底（Z=0），使用以下边界条件：
%   
%%  
%   $$ \frac{\partial \bar{P}}{\partial Z} = 0 \quad | Z=0^- $$
%   
%   $$ \frac{\partial y}{\partial Z} = 0 \quad | Z=0^- $$
%   
%   $$ \frac{\partial \bar{T}}{\partial Z} = 0 \quad | Z=0^- $$
%   
%%  
    P(1) = P(2) ;
    y(1) = y(2) ;
    T(1) = T(2) ;
%   
%%  
%   对于塔顶（Z=1），使用以下边界条件：
%   
%%  
%   $$ \frac{\partial \bar{P}}{\partial \tau} = \lambda(\bar{P_{eq}}-\bar{P}) \quad | Z=1^+ $$
%   
%   $$ y = y_{eq} \quad | Z=1^+ $$
%   
%   $$ \bar{T} = \bar{T}_{eq} \quad | Z=1^+ $$
%   
%%  
    y(N+2) = y_eq     ;
    T(N+2) = T_eq/T_0 ;
%   
%% 空间导数计算
%   对于一阶导数，使用加权基本无振荡（WENO）格式。
%   由于气体从塔顶流入塔内（速度为负），使用'downwind'格式。
%   
    Ph          = WENO(P, 'downwind')    ;
    dpdz(2:N+1) = (Ph(2:N+1)-Ph(1:N))/dz ;
    dpdzh(2:N)  = (P(3:N+1)-P(2:N))/dz   ;
    dpdzh(1)    =  2*(P(2)-P(1))/dz      ;
    dpdzh(N+1)  =  2*(P(N+2)-P(N+1))/dz  ;

    yh          = WENO(y, 'downwind')    ;
    dydz(2:N+1) = (yh(2:N+1)-yh(1:N))/dz ;
    
    Th          = WENO(T, 'downwind')    ;
    dTdz(2:N+1) = (Th(2:N+1)-Th(1:N))/dz ;
%   
%%  
%   *二阶导数*
%   
    d2ydz2(3:N) = (y(4:N+1)+y(2:N-1)-2*y(3:N))/dz/dz ;
    d2ydz2(2)   = (y(3)-y(2))/dz/dz                  ;
    d2ydz2(N+1) = (y(N)-y(N+1))/dz/dz                ;
    
    d2Tdz2(3:N) = (T(4:N+1)+T(2:N-1)-2*T(3:N))/dz/dz ;
    d2Tdz2(2)   =  4*(Th(2)+T(1)-2*T(2))/dz/dz       ;
    d2Tdz2(N+1) =  4*(Th(N)+T(N+2)-2*T(N+1))/dz/dz   ;
%   
%% 速度计算
%   根据Ergun方程计算气体速度。
%   
    ro_gh          = (P_0/R/T_0)*Ph(1:N+1)./Th(1:N+1)               ;
    viscous_term   = 150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2         ;
    kinetic_term_h = (ro_gh.*(MW_N2+(MW_CO2-MW_N2).*yh)).*(1.75*(1-epsilon)...
                     /2/r_p/epsilon)                                         ;
    vh = -sign(dpdzh).*(-viscous_term+(abs(viscous_term^2+4*kinetic_term_h... 
               .*abs(dpdzh)*P_0/L)).^(.5))/2./kinetic_term_h/v_0             ;
%   
%% 时间导数
%   
%   *1) 吸附相质量平衡*
%   
    q   = Isotherm(y, P*P_0, T*T_0, isotherm_params) ;
    q_1 = q(:, 1)*ro_s                               ;
    q_2 = q(:, 2)*ro_s                               ;
    
    k_1 = k_1_LDF*L/v_0 ;
    k_2 = k_2_LDF*L/v_0 ;
    
    dx1dt(2:N+1) = k_1*(q_1(2:N+1)/q_s0 - x1(2:N+1)) ;
    dx2dt(2:N+1) = k_2*(q_2(2:N+1)/q_s0 - x2(2:N+1)) ;
%   
%   *2) 塔的能量平衡*
%   
    sink_term = ((1-epsilon)*(ro_s*C_ps+q_s0*C_pa)+(epsilon.*ro_g(2:N+1).*C_pg)) ;
    
    % 2.1) 热传导
    transfer_term = K_z./v_0./L                             ;
    dTdt1(2:N+1)  = transfer_term.*d2Tdz2(2:N+1)./sink_term ;
    
    % 2.2) 对流
    PvT          = Ph(1:N+1).*vh(1:N+1)./Th(1:N+1) ;
    Pv           = Ph(1:N+1).*vh(1:N+1)            ;
    dTdt2(2:N+1) = -epsilon.*C_pg.*P_0./R./T_0.*((Pv(2:N+1)-Pv(1:N))- ... 
                    T(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz./sink_term      ;
    
    % 2.3) 吸附热
    generation_term_1 = (1-epsilon).*q_s0.*(-(deltaU_1-R*T(2:N+1)*T_0))./T_0 ;
    generation_term_2 = (1-epsilon).*q_s0.*(-(deltaU_2-R*T(2:N+1)*T_0))./T_0 ;
    dTdt3(2:N+1)      = (generation_term_1.*dx1dt(2:N+1)+...
                         generation_term_2.*dx2dt(2:N+1))./sink_term  ;
    
    % 2.4) 总和
    dTdt(2:N+1) = dTdt1(2:N+1) + dTdt2(2:N+1) + dTdt3(2:N+1) ;
%   
%   *3) 总质量平衡*
%   
    % 3.1) 对流
    dPdt1(2:N+1) = -T(2:N+1).*(PvT(2:N+1)-PvT(1:N))./dz  ;
    
    % 3.2) 吸附
    dPdt2(2:N+1) = -phi*T(2:N+1).*(dx1dt(2:N+1)+dx2dt(2:N+1)) ;
    
    % 3.3) 温度变化
    dPdt3(2:N+1) = P(2:N+1).*dTdt(2:N+1)./T(2:N+1)            ;
    
    % 3.4) 总和
    dPdt(2:N+1) = dPdt1(2:N+1) + dPdt2(2:N+1) + dPdt3(2:N+1)  ;
%   
%   *4) 组分质量平衡*
%   
    % 4.1) 扩散
    dydt1(2:N+1) = (1/Pe)*(d2ydz2(2:N+1)+(dydz(2:N+1).*dpdz(2:N+1)./P(2:N+1))... 
                  -(dydz(2:N+1).*dTdz(2:N+1)./T(2:N+1)))                        ;
    
    % 4.2) 对流
    ypvt         = yh(1:N+1).*Ph(1:N+1).*vh(1:N+1)./Th(1:N+1)        ;
    dydt2(2:N+1) = -(T(2:N+1)./P(2:N+1)).*((ypvt(2:N+1)-ypvt(1:N))... 
                   -y(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz             ;
    
    % 4.3) 吸附
    dydt3(2:N+1) = (phi*T(2:N+1)./P(2:N+1)).*((y(2:N+1)-1).*dx1dt(2:N+1)... 
                  + y(2:N+1).*dx2dt(2:N+1))                                ;
    
    % 4.4) 总和
    dydt(2:N+1) = dydt1(2:N+1) + dydt2(2:N+1) + dydt3(2:N+1) ;
%   
%%  边界导数
    dPdt(1)    = dPdt(2)                      ;
    dPdt(N+2)  = tau*L/v_0*(P_eq/P_0-P(N+2))  ;
    dydt(1)    = dydt(2)                      ;
    dydt(N+2)  = 0                            ; % Dirichlet BC for y
    dx1dt(1)   = 0                            ;
    dx2dt(1)   = 0                            ;
    dx1dt(N+2) = 0                            ;
    dx2dt(N+2) = 0                            ;
    dTdt(1)    = dTdt(2)                      ;
    dTdt(N+2)  = 0                            ; % Dirichlet BC for T
%   
%%  将导数导出到输出
    derivatives(1:N+2)        = dPdt(1:N+2)  ;
    derivatives(N+3:2*N+4)    = dydt(1:N+2)  ;
    derivatives(2*N+5:3*N+6)  = dx1dt(1:N+2) ;
    derivatives(3*N+7:4*N+8)  = dx2dt(1:N+2) ;
    derivatives(4*N+9:5*N+10) = dTdt(1:N+2)  ;
%   
end