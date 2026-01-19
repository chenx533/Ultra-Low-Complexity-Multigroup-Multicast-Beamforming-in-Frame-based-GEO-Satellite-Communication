function [minSE,Rate,SE] = getSumRate(prm,H,W,SNRdB)
% H：估计的信道状态信息矩阵，维度为[天线数×用户数×用户组数]
% W：Nt×B，天线数×波束数

[~,K,B] = size(H);
Rate = zeros(K,B);
for i=1:B
    temp = abs(W'*H(:,:,i)).^2;
    desire = temp(i,:);
    interf = sum(temp)-desire;
    Rate(:,i) = log2(1 + desire./(prm.xifit*(interf+1/prm.snr))).';
end
minSE = min(Rate);
%SE = K*sum(minSE);


% 不同信干噪比下的总频谱效率
snr = db2pow(SNRdB);                % 线性记法的信干噪比
SE = zeros(length(snr),1);

for j = 1:length(snr)
    for i=1:B
        temp = abs(W'*H(:,:,i)).^2;
        desire = temp(i,:);
        interf = sum(temp)-desire;
        SE(j) = SE(j)  + K*min(log2(1 + desire./(prm.xifit*(interf+1/snr(j)))));
    end
end

end