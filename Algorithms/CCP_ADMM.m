function W = CCP_ADMM(prm,H)
%CCP_ADMM算法：用于求解QoS问题
    
    % 参数提取和简化
    Nt = prm.nTxAntennas;               % 天线数
    B = prm.nGroup;                     % 多播组数
    Nu = prm.nUser;                     % 用户数
    sigma2 = (db2pow(prm.snr))^(-1);    % 噪声
    rho = 0.1;                          % ADMM惩罚参数
    Pn = 3;                             % 功率上限
    prm.xifit = 1.473;                  % 拟合因子
    th = 4*prm.xifit;                             % SINR阈值（QoS约束）
    Hres = reshape(H,[Nt,Nu*B]);        % 重新整合维度
    A = inv((2+rho)*eye(Nt) + rho*(Hres*Hres'));

    Wt = (1+1i)*eye(Nt,B)/sqrt(B);      % 初始化Wt
    W = zeros(Nt,B);
    
    tsca = 0;
    while(1)
    % for t = 1:prm.iternum  
        tsca = tsca + 1;
        %%%%%%%%%%%%%%%% 求解子问题Q(t) %%%%%%%%%%%%%%%%
        W0 = Wt;                            % 在子问题中初始化W
        Z = zeros(Nt,B);                    % 初始化z
        lambda = zeros(Nu*B,B);                         % 初始化λ
        Gamma = zeros(Nu*B,B);
    
        % 更新Gamma，按用户k子问题依次更新
        tadmm = 0;
        while (1)
            tadmm = tadmm + 1;
            for i = 1:Nu*B
                mk = ceil(i/Nu);
                hk = Hres(:,i);
                a = th*(sum(norm(hk'*W0-lambda(i,:))^2)-abs(hk'*W0(:,mk)-lambda(i,mk))^2);
                b = -2*abs(hk'*Wt(:,mk))^2;
                c = th*sigma2 - 2*real(Wt(:,mk)'*hk*(hk'*W0(:,mk)-lambda(i,mk))) + abs(hk'*Wt(:,mk))^2;
                if (a+c) < 0
                    pi = 0;
                else
                    p1 = b*th^2;
                    p2 = 2*b*th + c*th^2;
                    p3 = b + 2*c*th;
                    p4 = c + a;
                    rt = roots([p1 p2 p3 p4]);
                    id = ~imag(rt);
                    sre = rt(id);
                    pi = sre(sre>0);
                end
                Gamma(i,:) = (hk'*W0 - lambda(i,:))./(pi*th + 1);
                Gamma(i,mk) = pi*hk'*Wt(:,mk) + hk'*W0(:,mk) - lambda(i,mk);
            end
        
            % 更新v，按天线n子问题依次更新
            U = W0 - Z;
            for i = 1:Nt
                V(i,:) = min(sqrt(Pn)/norm(U(i,:)),1)*U(i,:);
            end
    
            % 更新w，直接用全局信息更新
            for i = 1:B
                W(:,i) = A*(rho*(V(:,i)+Z(:,i)) + rho*sum(Hres.*(Gamma(:,i).'+lambda(:,i).'),2));
            end
    
            % 更新对偶变量λ和z
            lambda = lambda + (Gamma-Hres'*W);
            Z = Z + (V-W);

            if (norm(W-W0)<0.001)||(tadmm>400)
                break
            else
                W0 = W;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (norm(W-Wt)<0.1)||(tsca>prm.iternum)   % CCP迭代停止条件
            break;
        else
            Wt = W;
        end
    end
end
