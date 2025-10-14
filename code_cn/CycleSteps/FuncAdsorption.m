function derivatives = FuncAdsorption(~, state_vars, params, isotherm_params)
%#codegen
%吸附：计算吸附步骤中状态变量的变化
%   与MATLAB ode求解器（ode15s）结合使用，求解变压吸附（PSA）步骤的偏微分方程组。
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
    
    y_0				=	params(23)	;
    P_inlet			=	params(26)	;
    ndot_0          =   P_0/R/T_0*v_0   ;
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
%%  
%   以下估算值为所使用的参数。轴向
%   弥散系数使用以下来自 Ruthvan "Pressure Swing Adsorption"
%   的方程计算
%   
%%  
%   $$ D_{L}= 0.7 D_{m} + r_{p} \cdot v_{0} $$
%   
%%  
%   气体浓度使用理想气体定律计算
%   
%%  
%   $$ C_{g}=\frac{\bar{P}}{R\bar{T}} \cdot \frac{P_{0}}{T_{0}} $$
%   
%% 计算所有使用的参数
    dz   = 1/N                                ;
    D_l  = 0.7*D_m + v_0*r_p                  ;
    Pe   = v_0*L/D_l                          ;
    phi  = R*T_0*q_s0*(1-epsilon)/epsilon/P_0 ;
    ro_g = P(1:N+2).*P_0/R./T(1:N+2)/T_0      ;
%   
%% 边界条件
%   吸附过程的边界条件，对应于塔两端
%   （轻组分和重组分产品端）均打开的情况。一股已知
%   温度和组成的原料气从重组分产品端送入塔内。
%   对于重组分产品端，使用以下边界条件：
%   
%%  
%   $$ \bar{P} = 1+\Delta \quad $$  或  $$ \quad \bar{v} = 1 \quad | Z=0^- $$
%   
%   $$ y = y_{0} \quad | Z=0^- $$
%   
%   $$ \bar{T} = \frac{T_{0}}{T_{0}} = 1 \quad | Z=0^- $$
%   
%%  
%   需要注意的是，对于入口压力，可以设置为一个常数，
%   或者将入口速度设置为常数，此时压力会发生变化以
%   满足入口处的速度约束。
    
    y(1) = y_0     ;
    T(1) = T_0/T_0 ;
    
    if params(end) == 1
        % 入口压力被指定为常数
        P(1) = P_inlet ;
    elseif params(end) == 0
        % 入口速度被指定为常数
        vh = zeros(N+1, 1) ;
        
        MW = MW_N2+(MW_CO2-MW_N2)*y(1) ;
    
        a_1   = 150*mu*(1-epsilon)^2*dz*L/2/4/r_p^2/epsilon^3/T(1)/T_0/R ;
        a_2_1 = 1.75*(1-epsilon)/2/r_p/epsilon/epsilon/epsilon*dz*L/2    ;
        a_2   = a_2_1/R/T(1)/T_0*ndot_0*MW                               ;
    
        a =  a_1+a_2             ;
        b =  P(2)/T(1)*P_0/R/T_0 ;
        c = -ndot_0              ;
    
        vh(1) = (-b+sqrt(b^2-4*a*c))/2/a/v_0 ;
    
        a_p = a_1*T(1)*T_0*R      ;
        b_p = a_2_1*MW/R/T(1)/T_0 ;
    
        P(1)  = ((a_p*vh(1)*v_0+P(2)*P_0)./(1-b_p*vh(1)*v_0*vh(1)*v_0))/P_0 ;
		
        % 为了重现Karson论文第五章的结果，
        % 必须添加下面的代码，这是我在Karson代码中发现的问题
        viscous_term=150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2;
        vh(1) = 1 ;
        P(1)=((viscous_term.*vh(1)*v_0*dz/2*L./P_0)+P(2))./(1-(dz/2*L/R/T(1)/T_0)*MW*(1.75*(1-epsilon)/2/r_p/epsilon)*vh(1).^2.*v_0.^2);
        
    else
        error('Please specify whether inlet velocity or pressure is constant for the feed step')
    end
%   
%%  
%   对于轻组分产品端，使用以下边界条件：
%   
%%  
%   $$ \bar{P} = 1 \quad | Z=1^+ $$
%   
%   $$ \frac{\partial y}{\partial \tau} = 0 \quad | Z=1^+ $$
%   
%   $$ \frac{\partial \bar{T}}{\partial \tau} = 0 \quad | Z=1^+ $$
%   
%%  
    y(N+2) = y(N+1)     ;
    T(N+2) = T(N+1)     ;
    if P(N+1) >= 1
        P(N+2) = 1      ;
    else
        P(N+2) = P(N+1) ;
    end
%   
%% 空间导数计算
%   
%   *一阶导数*
%   
%   对于一阶导数，每个体积单元的值是基于
%   其壁面上的值来估算的。这些值使用
%   加权基本无振荡（WENO）格式进行估算。
%   
%%  
%   $$ \frac{\partial f_j}{\partial Z}=\frac{f_{j+0.5}-f_{j-0.5}}{\Delta Z} $$
%   
%%  
%   压力：在体积单元中心和壁面上
    
    Ph          = WENO(P, 'upwind')      ;
    
    dpdz(2:N+1) = (Ph(2:N+1)-Ph(1:N))/dz ;
    dpdzh(2:N)  = (P(3:N+1)-P(2:N))/dz   ;
    dpdzh(1)    =  2*(P(2)-P(1))/dz      ;
    dpdzh(N+1)  =  2*(P(N+2)-P(N+1))/dz  ;
%   
%%  
%   摩尔分数：在体积单元中心
    
    yh          = WENO(y, 'upwind')      ;
    
    dydz(2:N+1) = (yh(2:N+1)-yh(1:N))/dz ;
%   
%%  
%   温度：在体积单元中心
    
    Th          = WENO(T, 'upwind')      ;
    
    dTdz(2:N+1) = (Th(2:N+1)-Th(1:N))/dz ;
%   
%%  
%   *二阶导数*
%   
%   二阶导数是基于节点上的值计算的。
%   只需要计算温度和摩尔分数的二阶导数。
%   需要注意的是，在塔的两端，没有发生扩散/传导，
%   所以在导数值中这些被设为0。
%   
%%  
%   $$ \frac{{\partial}^2 f_{j}}{\partial {Z}^2} = \frac{f_{j+1}+f_{j-1}-
%       2f_{j}}{{\Delta Z}^2} $$
%   
%%  
%   摩尔分数
    d2ydz2(3:N) = (y(4:N+1)+y(2:N-1)-2*y(3:N))/dz/dz ;
    d2ydz2(2)   = (y(3)-y(2))/dz/dz                  ;
    d2ydz2(N+1) = (y(N)-y(N+1))/dz/dz                ;
%   
%%  
%   温度
    d2Tdz2(3:N) = (T(4:N+1)+T(2:N-1)-2*T(3:N))/dz/dz ;
    d2Tdz2(2)   =  4*(Th(2)+T(1)-2*T(2))/dz/dz       ;
    d2Tdz2(N+1) =  4*(Th(N)+T(N+2)-2*T(N+1))/dz/dz   ;
%   
%% 速度计算
%   根据压力梯度计算气体在体积单元壁面处的间隙速度。
%   注意：请记住，Ergun方程中的速度是表观速度，
%   而不是间隙速度，这就是为什么在粘性项的分母中
%   是epsilon^2而不是epsilon^3，动能项也是如此，
%   其中是epsilon而不是epsilon^3。表观速度等于
%   间隙速度乘以空隙率。
%   
%%  
%   $$ U = v \cdot \varepsilon $$
%   
%   $$ - \frac{\partial \bar{P}}{\partial Z}\frac{P_0}{L} =
%   \frac{150\mu(1-\varepsilon)^2}{4r_{p}\varepsilon^2}\bar{v}v_0 +
%   \frac{1.75(1-\varepsilon)}{2r_{p}\varepsilon}
%   (\sum_{i}y_{i}MW_{i}C_{g})\bar{v}^2{v_{0}}^2 $$
%   
%%  
    ro_gh          = (P_0/R/T_0)*Ph(1:N+1)./Th(1:N+1)               ;
    
    viscous_term   = 150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2         ;
    kinetic_term_h = (ro_gh.*(MW_N2+(MW_CO2-MW_N2).*yh)).*(1.75*(1-epsilon)...
                     /2/r_p/epsilon)                                         ;
                    
    % 体积单元壁面上的速度
    vh = -sign(dpdzh).*(-viscous_term+(abs(viscous_term^2+4*kinetic_term_h... 
               .*abs(dpdzh)*P_0/L)).^(.5))/2./kinetic_term_h/v_0             ;
%   
%% 时间导数
%   
%%  
%   *1) 吸附相质量平衡（组分1和2的摩尔负载量）*
%       使用线性驱动力（LDF）模型计算气体的
%       吸附速率。假设LDF传质系数在
%       重要的压力和温度范围内是恒定的。
%   
%%  
%   $$ \frac{\partial x_{i}}{\partial \tau}q_{s0}= k_{i}({q_{i}}^*-q_{i}) $$
%   
%%  
%   1.1) 计算平衡摩尔负载量
    q   = Isotherm(y, P*P_0, T*T_0, isotherm_params) ;
    q_1 = q(:, 1)*ro_s                               ;
    q_2 = q(:, 2)*ro_s                               ;
%   
%%  
%   1.2) 计算LDF参数
    k_1 = k_1_LDF*L/v_0 ;
    k_2 = k_2_LDF*L/v_0 ;
%   
%%  
%   1.3) 计算时间导数
    dx1dt(2:N+1) = k_1*(q_1(2:N+1)/q_s0 - x1(2:N+1)) ;
    dx2dt(2:N+1) = k_2*(q_2(2:N+1)/q_s0 - x2(2:N+1)) ;
%   
%%  
%   *2) 塔的能量平衡（塔温）*
%       对于塔内的能量平衡，假设：
%       * 固相和气相之间处于热平衡状态
%       * 热传导同时通过固相和气相发生
%       * 对流只在气相中发生
%   
%%  
%   $$ \big[ \varepsilon C_{g}C_{p,g}+(1-\varepsilon)(C_{p,s}\rho_{s}+
%       C_{p,a}q_{s0})\big]\frac{\partial \bar{T}}{\partial \tau}= 
%       \frac{K_z}{v_0L} \frac{\partial^2 \bar{T}}{\partial Z^2}
%     - \varepsilon C_{g}C_{p,g} \bar{v}\frac{\partial \bar{T}}{\partial Z}
%     + \sum_{i} (1-\varepsilon)(-\Delta H_{i})\frac{q_{s0}}{T_{0}}\frac{
%       \partial \bar{x_{i}}}{\partial \tau} $$
%   
%%  
    % [J/m^3/K]
    sink_term = ((1-epsilon)*(ro_s*C_ps+q_s0*C_pa)+(epsilon.*ro_g(2:N+1).*C_pg)) ;
%   
%%  
%   2.1) 计算由于塔内热传导（固相和气相）
%        引起的温度变化
%   
%%  
%   $$ \frac{K_z}{v_0L} \frac{\partial^2 \bar{T}}{\partial Z^2} $$
%   
%%  
    transfer_term = K_z./v_0./L                             ;
    dTdt1(2:N+1)  = transfer_term.*d2Tdz2(2:N+1)./sink_term ;
%   
%%  
%   2.2) 计算由于对流引起的温度变化
%   
%%  
%   $$ -\varepsilon C_{g}C_{p,g} \bar{v}\frac{\partial \bar{T}}{\partial Z} $$
%   
%%  
    PvT          = Ph(1:N+1).*vh(1:N+1)./Th(1:N+1) ;
    Pv           = Ph(1:N+1).*vh(1:N+1)            ;
    dTdt2(2:N+1) = -epsilon.*C_pg.*P_0./R./T_0.*((Pv(2:N+1)-Pv(1:N))- ... 
                    T(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz./sink_term      ;
%   
%%  
%   2.3) 计算由于吸附/解吸焓变引起的温度变化
%   
%%  
%   $$ \sum_{i} (1-\varepsilon)(-\Delta
%   H_{i})\frac{q_{s0}}{T_{0}}\frac{\partial \bar{x_{i}}}{\partial \tau} $$
%   
%   $$ \Delta H_{i} = \Delta U_i - R\bar{T}T_{0} $$
%   
%%  
    generation_term_1 = (1-epsilon).*q_s0.*(-(deltaU_1-R*T(2:N+1)*T_0))./T_0 ;
    generation_term_2 = (1-epsilon).*q_s0.*(-(deltaU_2-R*T(2:N+1)*T_0))./T_0 ;
    
    dTdt3(2:N+1)      = (generation_term_1.*dx1dt(2:N+1)+...
                         generation_term_2.*dx2dt(2:N+1))./sink_term  ;
%   
%%  
%   2.4) 所有温度导数的总和
    dTdt(2:N+1) = dTdt1(2:N+1) + dTdt2(2:N+1) + dTdt3(2:N+1) ;
%   
%%  
%   *3) 总质量平衡*
%   
%%  
%   $$ \frac{\partial \bar{P}}{\partial \tau} = -\frac{\partial( \bar{v}
%      \bar{P}/\bar{T})}{\partial Z} - \Psi \bar{T} \sum_{i}\frac{\partial 
%      \bar{x_{i}}}{\partial \tau} + \frac{\bar{P}}{\bar{T}}\frac{\partial 
%      \bar{T}}{\partial \tau} $$
%   
%%  
%   3.1) 计算由于对流引起的压力变化
    dPdt1(2:N+1) = -T(2:N+1).*(PvT(2:N+1)-PvT(1:N))./dz  ;
%   
%%  
%   3.2) 计算由于吸附/解吸引起的压力变化
    dPdt2(2:N+1) = -phi*T(2:N+1).*(dx1dt(2:N+1)+dx2dt(2:N+1)) ;
%   
%%  
%   3.3) 计算由于温度变化引起的压力变化
    dPdt3(2:N+1) = P(2:N+1).*dTdt(2:N+1)./T(2:N+1)            ;
%   
%%  
%   3.4) 所有压力变化的总和
    dPdt(2:N+1) = dPdt1(2:N+1) + dPdt2(2:N+1) + dPdt3(2:N+1)  ;
%   
%%  
%   *4) 组分质量平衡（基于摩尔分数）*
%   
%%  
%   $$ \frac{\partial y}{\partial \tau} = \frac{1}{Pe} \big(\frac{{
%      \partial}^2 y}{\partial {Z}^2}+\frac{1}{\bar{P}}\frac{\partial
%      \bar{P}}{\partial Z}\frac{\partial y}{\partial Z}-\frac{1}{\bar{T}}
%      \frac{\partial \bar{T}}{\partial Z}\frac{\partial y}{\partial Z}
%      \big)-\bar{v}\frac{\partial y}{\partial Z}+\frac{\Psi \bar{T}}{
%      \bar{P}} \big((y-1)\frac{\partial \bar{x_{1}}}{\partial \tau}+y
%      \frac{\partial\bar{x_{2}}}{\partial \tau}\big) $$
%   
%%  
%   4.1) 计算由于扩散引起的摩尔分数变化
%   
%%  
%   $$ \frac{1}{Pe} \big(\frac{{\partial}^2 y}{\partial {Z}^2}+\frac{1}{
%      \bar{P}}\frac{\partial \bar{P}}{\partial Z}\frac{\partial y}{
%      \partial Z}-\frac{1}{\bar{T}}\frac{\partial \bar{T}}{\partial Z}
%      \frac{\partial y}{\partial Z} \big) $$
%   
%%  
    dydt1(2:N+1) = (1/Pe)*(d2ydz2(2:N+1)+(dydz(2:N+1).*dpdz(2:N+1)./P(2:N+1))... 
                  -(dydz(2:N+1).*dTdz(2:N+1)./T(2:N+1)))                        ;
%   
%%  
%   4.2) 计算由于对流引起的摩尔分数变化
%   
%%  
%   $$ -\bar{v}\frac{\partial y}{\partial Z} $$
%   
%%  
    ypvt         = yh(1:N+1).*Ph(1:N+1).*vh(1:N+1)./Th(1:N+1)        ;
    dydt2(2:N+1) = -(T(2:N+1)./P(2:N+1)).*((ypvt(2:N+1)-ypvt(1:N))... 
                   -y(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz             ;
%   
%%  
%   4.3) 计算由于吸附/解吸引起的摩尔分数变化
%   
%%  
% $$ \frac{\Psi \bar{T}}{\bar{P}} \big((y-1)\frac{\partial \bar{x_{1}}}{
%    \partial \tau}+y\frac{\partial \bar{x_{2}}}{\partial \tau}\big) $$
%   
%%  
    dydt3(2:N+1) = (phi*T(2:N+1)./P(2:N+1)).*((y(2:N+1)-1).*dx1dt(2:N+1)... 
                  + y(2:N+1).*dx2dt(2:N+1))                                ;
%   
%%  
%   4.4) 所有摩尔分数变化的总和
    dydt(2:N+1) = dydt1(2:N+1) + dydt2(2:N+1) + dydt3(2:N+1) ;
%   
%%  边界导数
    dPdt(1)    = 0         ;
    dPdt(N+2)  = 0         ;
    dydt(1)    = 0         ;
    dydt(N+2)  = dydt(N+1) ;
    dx1dt(1)   = 0         ;
    dx2dt(1)   = 0         ;
    dx1dt(N+2) = 0         ;
    dx2dt(N+2) = 0         ;
    dTdt(1)    = 0         ;
    dTdt(N+2)  = dTdt(N+1) ;
%   
%%  将导数导出到输出
    derivatives(1:N+2)        = dPdt(1:N+2)  ;
    derivatives(N+3:2*N+4)    = dydt(1:N+2)  ;
    derivatives(2*N+5:3*N+6)  = dx1dt(1:N+2) ;
    derivatives(3*N+7:4*N+8)  = dx2dt(1:N+2) ;
    derivatives(4*N+9:5*N+10) = dTdt(1:N+2)  ;
%   
end