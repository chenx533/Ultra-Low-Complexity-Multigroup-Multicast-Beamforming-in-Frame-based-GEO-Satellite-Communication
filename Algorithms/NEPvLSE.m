function W = NEPvLSE(prm,H,PAalg)
    arguments
        prm struct
        H (:,:,:) double
        PAalg {mustBeMember(PAalg,["LDD","WF","none"])} = "none";
    end

    Nt = prm.nTxAntennas;
    B = prm.nGroup;
    Nu = prm.nUser;
    tau = 100;

    IB = eye(B);
    A = zeros(Nt*B,Nt*B,Nu,B);
    C = zeros(Nt*B,Nt*B,Nu,B);
    for b = 1:B
        for m = 1:Nu
            h = H(:,m,b);
            G = h*h';
            A(:,:,m,b) = kron(IB,G) - (prm.xifit-1)*kron(diag(IB(:,b)),G)/prm.xifit + db2pow(prm.snr)^(-1)*eye(B*Nt);
            C(:,:,m,b) = kron(IB,G) - kron(diag(IB(:,b)),G) + db2pow(prm.snr)^(-1)*eye(B*Nt);
        end
    end
    % 创建非线性函数的句柄
    
    R = @(w) log2(pagemtimes(pagemtimes(w',A),w)./pagemtimes(pagemtimes(w',C),w));
    %f = @(w) real(reshape(-log(sum((pagemtimes(pagemtimes(w',A),w)./pagemtimes(pagemtimes(w',C),w)).^(-tau/log(2)),3))/tau,[B,1]));
    coff = @(w) exp(-tau*R(w))./repmat(sum(exp(-tau*R(w)),3),[1,1,Nu,1]);
    
    M = @(w) inv(sum(coff(w).*(C./pagemtimes(pagemtimes(w',C),w)),[3,4])) * sum(coff(w).*(A./pagemtimes(pagemtimes(w',A),w)),[3,4]); % 不动点迭代法
    
    % 迭代求解
    w = randn(Nt*B,1,"like",1i);
    tnepv = 0;
    %while (1)
    for t = 1:prm.iternum
        w_1 = w;
        w = M(w)*w;
        %w = w/norm(w);
        tnepv = tnepv + 1;
        if norm(w-w_1)<0.1    % 迭代变化情况
            break
        end
    end
    w = w/norm(w);
    W = reshape(w,[Nt,B]);

    % 功率分配
    SE = reshape(R(w),[Nu,B]);
    [~,index]=min(SE);
    
    Gmin = zeros(Nt,B);
    for i = 1:B
        Gmin(:,i) = H(:,index(i),i);
    end
    G = abs(Gmin'*W).^2;
    switch PAalg
        case "LDD"
            p = powerAllocate_Lagrange(prm,G);
        case "WF"
            p = powerAllocate_WaterFilling(prm,G);
        case "none"
            return
    end
    W = W.*sqrt(p).';   % 总的波束成形矩阵
end



function p = powerAllocate_Lagrange(prm,H)
    B = prm.nGroup;             % 波束数量
    p = 2*ones(B,1);              % 功率分配系数
    Pt = sum(p);
    lambda = ones(B,1);         % 初始拉格朗日乘子
    hbb = diag(H);              % 期望用户的信道
    sigma2 = (db2pow(prm.snr))^(-1);    % 噪声
    th = 4;                   % SINR阈值

    options = optimset('FunValCheck','off');
    
    gamma = (p.*hbb)./(prm.xifit*(H*p-p.*hbb+sigma2));  % 初始值gamma_0
    y = sqrt((gamma+1).*hbb.*p)./(hbb.*p+prm.xifit*(H*p-p.*hbb+sigma2));    % 初始值y_0
    
    tldd = 0;      % 迭代次数下标
    while (1)
        ypow = y.^2;            % 记录y^2的临时变量
        for i = 1:B
            if hbb(i)*p(i)-th*prm.xifit*(H(i,:)*p - hbb(i)*p(i) + sigma2) > 0.01
                lambda(i) = 0;
            else
                hjm = H(:,i);
                hjm(i) = 0;
                lambda(i) = ypow(i)-(-sqrt((hbb(i)^2)*ypow(i)*(1+gamma(i))./(th*prm.xifit*(H(i,:)*p - hbb(i)*p(i) + sigma2)))-prm.xifit*(ypow + th*lambda).'*hjm)/hbb(i);
            end
        end

        Hjm = H-diag(diag(H));
        p = (ypow.*(1+gamma).*hbb)./((ypow-lambda).*hbb+prm.xifit*((ypow + th*lambda).'*Hjm).').^2;
        
        if(sum(p)>Pt)
            % 归一化功率
            % index = find(lambda==0);
            % temp=(B-sum(p(lambda~=0)))*p(index)/sum(p(index));
            % p(index) = temp;
            weight = (lambda~=0);
            Pmax = Pt-weight.'*p;
            if sum(~weight)~=0
                findv = @(v) (~weight.')*((ypow.*(1+gamma).*hbb)./((ypow-lambda).*hbb+prm.xifit*((ypow + th*lambda).'*Hjm).'+v).^2)-Pmax;
                v = fzero(findv,0,options);
            end
            temp = (ypow.*(1+gamma).*hbb)./((ypow-lambda).*hbb+prm.xifit*((ypow + th*lambda).'*Hjm).'+v).^2;
            p(lambda==0)=temp(lambda==0);
        end

        gamma1 = (p.*hbb)./(prm.xifit*(H*p-p.*hbb+sigma2));

        y1 = sqrt((gamma+1).*hbb.*p)./(hbb.*p+prm.xifit*(H*p-p.*hbb+sigma2));
        if (norm(y1-y) < 0.001 || tldd>100)
            break;
        else
            y = y1;
            gamma = gamma1;
        end
        tldd = tldd + 1;
    end
end



function p = powerAllocate_WaterFilling(prm,H)
    B = prm.nGroup;
    sigma2 = (db2pow(prm.snr))^(-1);    % 噪声
    Pt = 14;
    th = 4;
    [hbb,ind] = sort(diag(H),'descend');
    lambda = prm.xifit*sigma2./hbb;
    tau = lambda*th;
    for i = 1:B
        index = B-i+1;
        nu = (B-i+1)/(Pt - sum(tau(index+1:end)) + sum(lambda(1:index)));
        p(1:index) = 1./nu-lambda(1:index);
        if all(p(1:index)>=tau(1:index))
            break
        end
        if p(index)<tau(index)
            p(index:end)=tau(index:end);
        end
    end
    p = p(ind);
end