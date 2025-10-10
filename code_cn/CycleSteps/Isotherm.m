function q=Isotherm(y, P, T, isotherm_parameters)
%hypotheticalisotherm- 计算双组分混合物的摩尔负载量
%   计算双组分混合物的总负载量。
%   假设该等温线能用双活性位Langmuir竞争性吸附等温线表示，
%   且其参数与温度相关。
%   该等温线可以基于浓度 [mol/m^3] 或分压 [Pa] 表示。
%   输入:
%   y: 组分一的摩尔分数 [-]。假设组分二的摩尔分数为1-y
%   P: 气体的总压力 [Pa]
%   T: 气体的温度 [K]
%   isothermparams: 等温线参数。其结构请参见 Input_PSA 脚本
%   input_units: 指定等温线是基于浓度还是分压

    R=8.314;

    q_s_b_1=isotherm_parameters(1);
    q_s_d_1=isotherm_parameters(3);
    q_s_b_2=isotherm_parameters(2);
    q_s_d_2=isotherm_parameters(4);
    b_1=isotherm_parameters(5);
    d_1=isotherm_parameters(7);
    b_2=isotherm_parameters(6);
    d_2=isotherm_parameters(8);
    deltaU_b_1=isotherm_parameters(9);
    deltaU_d_1=isotherm_parameters(11);
    deltaU_b_2=isotherm_parameters(10);
    deltaU_d_2=isotherm_parameters(12);


    B_1=b_1*exp(-deltaU_b_1/R./T);
    D_1=d_1*exp(-deltaU_d_1/R./T);
    B_2=b_2*exp(-deltaU_b_2/R./T);
    D_2=d_2*exp(-deltaU_d_2/R./T);

    if isotherm_parameters(13) == 0

        P_1=y.*P;
        P_2=(1-y).*P;
        input_1=P_1;
        input_2=P_2;

    elseif isotherm_parameters(13) ==1
        C_1=y.*P./R./T;
        C_2=(1-y).*P./R./T;
        input_1=C_1;
        input_2=C_2;
    else
        error('Please specify whether the isotherms are in terms of Concentration or Partial Pressure')

    end


    q1_b=q_s_b_1.*B_1.*input_1./(1+B_1.*input_1+B_2.*input_2);
    q1_d=q_s_d_1.*D_1.*input_1./(1+D_1.*input_1+D_2.*input_2);

    q1=q1_b+q1_d;

    q2_b=q_s_b_2.*B_2.*input_2./(1+B_1.*input_1+B_2.*input_2);
    q2_d=q_s_d_2.*D_2.*input_2./(1+D_1.*input_1+D_2.*input_2);

    q2=q2_b+q2_d;
	
	q = [q1, q2] ;

end