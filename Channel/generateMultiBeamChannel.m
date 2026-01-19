%% generate multi-beam channel
% 用途：根据卫星多播场景布局和信道参数，生成多波束卫星信道

% 输入：prm —— 系统参数及卫星参数，形式为结构体

% 输出：H —— 估计的信道状态信息矩阵，不带相位噪声，维度为[反射面天线数×用户数×用户组数]

% 执行示例：
% H = generateMultiBeamChannel(prm,user);

% 最后修改日期：2024-06-03
% =========================================================================
function H = generateMultiBeamChannel(prm,user)
    % 计算大尺度参数
    psi = (prm.lightSpeed/(4*pi*prm.carrierFrequency*prm.satHeight))^2/( ...
        prm.boltzmann*prm.bandWidth);

    % 远场波束夹角
    userAzEl = zeros(2,prm.nUser,prm.nGroup);
    groupAzEl = zeros(2,prm.nGroup);                            
    for i = 1:prm.nGroup
        % 将群组方位角和俯仰角转换至相控阵视角
        groupAzEl(1,i) = 90 - prm.groupAngle(1,i);
        groupAzEl(2,i) = mod(prm.groupAngle(2,i)-180,360);
        % 将用户方位角和俯仰角转换至相控阵视角
        for n = 1:prm.nUser
            userAzEl(1,n,i) = 90 - user{i,n}.elevationAngle;
            userAzEl(2,n,i) = mod(user{i,n}.azimuthAngle-180,360);
            % [userAzEl(1,n,i),userAzEl(2,n,i)]=glb2loc(user{i,n}.azimuthAngle ...
            %     -180,-user{i,n}.elevationAngle);
        end
    end
    varPhi = getOffAxisAngle(groupAzEl,userAzEl);   % 计算用户与各波束夹角
    U = 2.07123*sind(varPhi)/sind(prm.angle3db);
    Gr = db2pow(prm.rxGT)*ones(prm.nGroup,1);
    Gt = db2pow(prm.txGain-30)*ones(prm.nGroup,1);  
    zeroind = (U==0);                               % 角度为0的地方的方向图索引
    pattern = (besselj(1,U)./(2*U) + 36*besselj(3,U)./(U.^3)).^2;
    pattern(zeroind) = 1;                           % 角度为0的地方方向图增益为1
    B = repmat(Gr.*Gt,1,prm.nUser,prm.nGroup).*pattern; % 计算期望波束和干扰波束增益
    
    % 雨衰系数
    rainFadingdb = -2.6+sqrt(1.63)*randn(prm.nTxAntennas,prm.nUser, ...
        prm.nGroup);              % 雨衰服从对数正态分布
    r = db2pow(rainFadingdb);     % 将雨衰从对数记法转换为线性记法

    % 信道相位向量，生成[0,2pi]范围内的随机相位
    theta = 2*pi*rand(prm.nTxAntennas,prm.nUser,prm.nGroup);

    % 总信道矩阵根据大尺度衰落，波束增益，雨衰和随机相位生成
    H = sqrt(psi).*sqrt(B).*sqrt(r).*exp(1i*theta);
end