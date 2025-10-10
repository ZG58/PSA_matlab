function flux_w = WENO(flux_c, FlowDir)
%WENO: 应用加权基本无振荡格式
%   接收有限体积单元中心处的通量和流动方向（顺风 -
%   逆风），流动方向由速度或压降决定。
%   如果 v>=0 或 dP<0，则应用顺风格式。如果 v<0 或 dP>0，则
%   应用逆风格式。该函数返回有限体积
%   单元壁面处的通量。顺风（upwind）意味着计算在域上从左到右进行，
%   或者就PSA塔而言，是"相对于进料口的并流流动"。
%   逆风（downwind）意味着计算在域上从右到左进行，
%   或者就PSA塔而言，是"相对于进料口的逆流流动"。
%   
%   输入:
%       flux_c : 有限体积单元中心处的通量
%       FlowDir: 流动方向。选项: upwind (顺风) 和 downwind (逆风)
%   
%   输出:
%       flux_w : 有限体积单元边缘或壁面处的通量
%   
%%  
%   对于并流（顺风），有限体积单元壁面处的通量
%   计算如下：
%   
%%  
%   $$ f_{j+0.5}=\frac{\alpha_{0,j}}{\alpha_{0,j}+\alpha_{1,j}} 
%   \Big[\frac{1}{2}(f_{j}+f_{j+1})\Big] + \frac{\alpha_{1,j}}{\alpha_{0,j}+
%   \alpha_{1,j}}\Big[\frac{3}{2}f_{j}-\frac{1}{2}f_{j-1}\Big] $$
%   
%   $$ \alpha_{0,j}= \frac{2/3}{(f_{j+1}-f_{j}+\delta)^4} $$
%   
%   $$ \alpha_{1,j}= \frac{1/3}{(f_{j}-f_{j-1}+\delta)^4} $$
%   
%%  
%   对于逆流（逆风），有限体积单元壁面处的通量
%   计算如下：
%   
%%  
%   $$ f_{j+0.5}=\frac{\alpha_{0,j}}{\alpha_{0,j}+\alpha_{1,j}} 
%   \Big[\frac{1}{2}(f_{j}+f_{j+1})\Big] + \frac{\alpha_{1,j}}{\alpha_{0,j}+
%   \alpha_{1,j}}\Big[\frac{3}{2}f_{j+1}-\frac{1}{2}f_{j+2}\Big] $$
%   
%   $$ \alpha_{0,j}= \frac{2/3}{(f_{j}-f_{j+1}+\delta)^4} $$
%   
%   $$ \alpha_{1,j}= \frac{1/3}{(f_{j+1}-f_{j+2}+\delta)^4} $$
%   
%%  
    oo     = 10^-10              ;
    [N, m] = size(flux_c)        ;
    N      = N-2                 ;
    flux_w = zeros(N+1, m)       ;
    alpha0 = zeros(size(flux_c)) ;
    alpha1 = zeros(size(flux_c)) ;
    
    % 域边界处的通量
    flux_w(1, :)   = flux_c(1, :)   ;
    flux_w(N+1, :) = flux_c(N+2, :) ;
    
    if strcmpi(FlowDir, 'upwind') == 1
        
        alpha0(2:N, :) =(2/3)./((flux_c(3:N+1, :)-flux_c(2:N, :)+oo).^4) ;
        alpha1(3:N, :) =(1/3)./((flux_c(3:N, :)-flux_c(2:N-1, :)+oo).^4) ;
        alpha1(2, :)   =(1/3)./((2*(flux_c(2, :)-flux_c(1, :))+oo).^4)   ;
        
        flux_w(3:N, :) = (alpha0(3:N, :)./(alpha0(3:N, :)+alpha1(3:N, :)))...      
                       .*((flux_c(3:N, :)+flux_c(4:N+1, :))./2)+(alpha1(3:N, :)... 
                       ./(alpha0(3:N, :)+alpha1(3:N, :))).*(1.5*flux_c(3:N, :)...  
                        -.5*flux_c(2:N-1, :))                                     ;
        
        flux_w(2, :)   = (alpha0(2, :)./(alpha0(2, :)+alpha1(2, :)))...    
                       .*((flux_c(2, :)+flux_c(3, :))./2)+(alpha1(2, :)... 
                       ./(alpha0(2, :)+alpha1(2, :))).*(2*flux_c(2, :)...  
                         -flux_c(1, :))                                   ;
        
    elseif strcmpi(FlowDir, 'downwind') == 1
        
        alpha0(2:N, :)   = (2/3)./((flux_c(2:N, :)-flux_c(3:N+1, :)+oo).^4)   ;
        alpha1(2:N-1, :) = (1/3)./((flux_c(3:N, :)-flux_c(4:N+1, :)+oo).^4)   ;
        alpha1(N, :)     = (1/3)./((2*(flux_c(N+1, :)-flux_c(N+2, :))+oo).^4) ;
        
        flux_w(2:N-1, :) = (alpha0(2:N-1, :)./(alpha0(2:N-1, :)+alpha1(2:N-1, :)))...   
                         .*((flux_c(2:N-1, :)+flux_c(3:N, :))./2)+(alpha1(2:N-1, :)...  
                         ./(alpha0(2:N-1, :)+alpha1(2:N-1, :))).*(1.5*flux_c(3:N, :)... 
                          -.5*flux_c(4:N+1, :))                                        ;
                      
        flux_w(N, :)     = (alpha0(N, :)./(alpha0(N, :)+alpha1(N, :)))...      
                         .*((flux_c(N, :)+flux_c(N+1, :))./2)+(alpha1(N, :)... 
                         ./(alpha0(N, :)+alpha1(N, :))).*(2*flux_c(N+1, :)...  
                           -flux_c(N+2, :))                                   ;
    else
        error('Please specify the direction of flow. OPTIONS: upwind and downwind')
    end 
%   
end