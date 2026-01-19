function w = SDRHeuristic(prm,H)
    Nt = prm.nTxAntennas;
    K = prm.nUser;
    M = prm.nGroup;
    Pt = 2;
    delta = 0.1;
    sigma2 = (db2pow(prm.snr))^(-1);    % 噪声
    X = getChannelCorr(prm,H);    % 获取二阶信道矩阵，维度为[天线数×天线数×用户数×用户组数]
    [H_qos_des, H_qos_int] = getQoSChannMat(X);     % 生成便于CVX计算的信道矩阵
    th = 4;

    w = 1*(1+1i)*eye(Nt,M);
    p0 = zeros(1,M);
    while(1)
    % for t = 1:prm.iternum 
        %-----------------------------计算γ--------------------------------
        [minSE,Rate,~] = getSumRate(prm,H,w,prm.snr);
        gamma = repelem(2.^(minSE)-1,K).';
        %-----------------------------求解w--------------------------------
        cvx_clear
        cvx_solver SeDuMi
        cvx_begin sdp quiet
            variable W(Nt,Nt,M) hermitian
            minimize trace(sum(W,3))
            subject to
                real(H_qos_des*vec(W)) - prm.xifit*max(gamma,th).*(real(H_qos_int*vec(W)-H_qos_des*vec(W)) + sigma2) >= zeros(K*M,1);
                W >= 0;
        cvx_end
        w = recover_weight(W);
        %-----------------------------计算p--------------------------------
        p = vecnorm(w).^2;
        s = log(p);
        v = w./sqrt(p);
        %-----------------------------选择k--------------------------------
        [~,index] = min(Rate);
        Gmin = zeros(Nt,M);
        for i = 1:M
            Gmin(:,i) = H(:,index(i),i);
        end
        V = abs(Gmin'*v).^2;
        %-----------------------------计算r--------------------------------
        z = sigma2./diag(V);
        nu = V./diag(V) - eye(M);
        I = sum(p.*nu,2) + z;
        sinr = p.'./I;
        g = -1./(sinr.*I);
        r = p.'.*(g - sum(nu.*sinr.*g,2));
        s = s - (delta*r).';
        p = exp(s);
        p = p/sum(p)*Pt;
        %-----------------------------更新w--------------------------------
        w = sqrt(p).*v;
        if sum(abs(p-p0))<0.1
            break
        else
            p0 = p;
        end
    end
    %w = w/norm(w,'fro')*sqrt(P);
end

function X = getChannelCorr(prm,H)
    Nt = prm.nTxAntennas;
    K = prm.nUser;
    M = prm.nGroup;
    X = zeros(Nt,Nt,K,M);
    for i = 1:M
        for u = 1:K
            h = H(:,u,i);
            corr = h*h';
            X(:,:,u,i) = corr;
        end
    end
end

function [H_qos_des, H_qos_int] = getQoSChannMat(X)
    [Nt,~,K,M] = size(X);
    H_qos_des = [];
    H_qos_int = zeros(K*M, Nt*Nt*M);
    for i = 1:M
        temp_des = reshape(X(:,:,:,i), [Nt*Nt,K]).';
        temp_int = repmat(temp_des,[1,M]);
        H_qos_des = blkdiag(H_qos_des,temp_des);
        index = (i-1)*K+1:i*K;
        H_qos_int(index,:) = temp_int;
    end
end

function w = recover_weight(W)
    [Nt,~,M]=size(W);
    w=zeros(Nt,M);
    for i=1:M
        [V,D]=eig(W(:,:,i));
        [d,ind] = sort(diag(D));
        U=V(:,ind);
        w(:,i)=sqrt(d(end))*U(:,end);
    end
end