function W = MGMBF(prm,H)
% 基于CCCP的多播波束成形算法

    % 参数提取和简化
    Nt = prm.nTxAntennas;               % 天线数
    B = prm.nGroup;                     % 多播组数
    Nu = prm.nUser;                     % 用户数
    sigma2 = (db2pow(prm.snr))^(-1);    % 噪声
    th = 4;                             % SINR阈值（QoS约束）
    Pt = 2;                             % 功率上限
    prm.xifit = 1.473;                  % 拟合因子

    % 辅助矩阵
    extract = ones(Nu*B,B) - kron(eye(B),ones(Nu,1));

    % CVX求解过程用到的变形信道矩阵
    H_2D = reshape(H,[Nt,Nu*B]);
    H_rsp1 = H_2D';
    H_rsp2 = permute(H_2D,[1,3,2]);

    % 初始化W和gamma
    W0 = (1+1i)*eye(Nt,B)/sqrt(B);    % 波束成形矩阵初始值
    usrSINR = (sum_square_abs(H_rsp1*W0,2) + sigma2)./(sum_square_abs((H_rsp1*W0).*extract,2) + sigma2) - 1;         % gamma初始值
    sinr0 = min(reshape(usrSINR,[Nu,B])).';

    % SCA迭代
    while(1)
    % for t = 1:prm.iternum    
        H_rsp1W0 = H_rsp1*W0;           % 辅助计算的中间变量
        sinr0s = repelem(sinr0,Nu);     % 按用户数扩展的gamma

        % 更新展开点的函数值F(t-1)
        F = (sum_square_abs(H_rsp1W0,2) + sigma2)./(1 + sinr0s/prm.xifit);
        
        % 计算nablaF(t-1)
        a = permute(2*H_rsp1W0./(1 + sinr0s/prm.xifit),[3,2,1]);
        nFW = reshape(pagemtimes(H_rsp2,a),[Nt*B,Nu*B]);        % 计算对W的导数部分
        %------------------------------------------------
        temp = reshape(-F./(1 + sinr0s/prm.xifit),[Nu,B]).';    % 计算对gamma的导数部分
        nFgam = [];
        for i = 1:B  nFgam=blkdiag(nFgam,temp(i,:));  end
        %------------------------------------------------
        nablaF = [nFW;nFgam];   % 更新一阶导数矩阵

        % CVX求解
        cvx_clear;
        cvx_solver SeDuMi
        cvx_expert true;
        cvx_begin quiet
            variable W(Nt,B) complex
            variable sinr(B,1)
            maximize Nu*sum_log(1+sinr/prm.xifit)/log(2)
            subject to
                % P4约束C9
                (sum_square_abs((H_rsp1*W).*extract,2) + sigma2) - F - real(nablaF'*[vec(W)-vec(W0);sinr-sinr0]) <= 0;  % F为展开点的函数值，nablaF为F矩阵函数对变量W和y的一阶导
                % P3约束C8
                sinr >= th*prm.xifit;
                % P1约束C5
                sum_square_abs(vec(W)) <= Pt;
        cvx_end

        if abs(sum_log(1+sinr/prm.xifit)-sum_log(1+sinr0/prm.xifit))<=0.1
            return
        else
            W0 = W;     % 更新W的展开点
            sinr0 = sinr;   % 更新y的展开点
        end
    end
end