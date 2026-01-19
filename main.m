clear all;
clc;

seed = randi(500);  
disp(seed);
rng(seed);    %seed=262(2个),106(0个)复现论文中的最小SE柱状图225,395

% 场景参数 
prm.lightSpeed = physconst('LightSpeed');             % 光速
prm.satHeight = 3.5768e7;                             % 卫星轨道高度（km）
prm.nGroup = 7;                                       % 多播组数量
prm.nUser = 5;                                        % 每组用户数量
prm.groupAngle = [90 88.23 88.23 88.23 88.23 88.23 88.23;0 0 60 120 180 240 300]; % 波束离轴角按卫星小区形状排布，见programming idea
user = cell(prm.nGroup,prm.nUser);                    % 用于存放用户参数的结构体数组
for i = 1:prm.nGroup                                  % 生成场景布局
    for n = 1:prm.nUser
        user{i,n}.elevationAngle = prm.groupAngle(1,i) + 0.05*rand; % 设置用户俯仰角
        user{i,n}.azimuthAngle = prm.groupAngle(2,i) + 0.05*rand;   % 设置用户方位角
    end
end
        
% 星地信道参数
prm.bandWidth = 20e6;                                 % 设置带宽，默认为50MHz
prm.carrierFrequency = 2e9;                           % 载波频率为2GHz，S波段
prm.lambda = prm.lightSpeed/prm.carrierFrequency;     % 计算波长
prm.phaseNoiseVar = 10*pi/180;                        % 信道相位噪声偏差为10°
prm.rainFadingMean = -2.6;                            % 雨衰服从对数分布，均值为-2.6dB
prm.rainFadingVar = 1.63;                             % 雨衰服对数分布方差为1.63dB
prm.boltzmann = physconst("Boltzmann");               % 玻尔兹曼常数为1.38×10^-23 J/K

% 用户端接收机参数
prm.noise = 1;                                        % 接收机归一化噪声系数
prm.rxGT = 19;                                        % 接收机的G/T值（dB/K）

% 卫星发射机参数
prm.nTxAntennas = prm.nGroup;                         % 设置卫星发射天线数与波束数相同
prm.txGain = 54;                                      % 发射天线增益（dBi）
prm.angle3db = 0.4011;                                % 卫星发射天线的3dB波束角
prm.snr = 10;                                         % 信噪比
prm.xifit = 1.473;                                    % DVB-S2X拟合系数
prm.iternum = 15;

%------------------------------卫星多波束信道---------------------------------
H = generateMultiBeamChannel(prm,user);
%H = load("channMatrix_Monte_1.mat","Hgg");
%-------------------------------波束成形算法----------------------------------
try
    tic;W1 = NEPvLSE(prm,H,'LDD');toc;
    W2 = NEPvLSE(prm,H,'WF');
    W3 = MGMBF(prm,H);
    W4 = SDRHeuristic(prm,H);
    tic;W = CCP_ADMM(prm,H);toc;
    % --------------------------------性能评估------------------------------------
    [minSE1,Rate1,SE1] = getSumRate(prm,H,W1,prm.snr);
    [minSE2,Rate2,SE2] = getSumRate(prm,H,W2,prm.snr);
    [minSE3,Rate3,SE3] = getSumRate(prm,H,W3,prm.snr);
    [minSE4,Rate4,SE4] = getSumRate(prm,H,W4,prm.snr);
    [minSE,Rate,SE] = getSumRate(prm,H,W,prm.snr);
catch
    disp('原问题不可行或不精确解');
end
