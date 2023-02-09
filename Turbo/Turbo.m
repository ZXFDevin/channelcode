%Turbo码仿真测试
clear
clc

SNR= [-6, -5.75, -5.5, -5.25];                                  %设定信噪比
frameLoop = [1*(10^4),2*(10^4),5*(10^4),1*(10^5)];                %每个信噪比处仿真的帧数

frameErrorRate = zeros(1,length(SNR));  


frmLen = 256;                                   %设定每一帧信息数量
rng default

intrlvrIndices = randperm(frmLen);

hTEnc = comm.TurboEncoder('TrellisStructure',poly2trellis(4, ...
    [13 15 17],13),'InterleaverIndices',intrlvrIndices);

hTDec = comm.TurboDecoder('TrellisStructure',poly2trellis(4, ...
    [13 15 17],13),'InterleaverIndices',intrlvrIndices, ...
    'NumIterations',4);

hMod = comm.BPSKModulator;



for snr_index = 1:1:length(SNR)
    %打印中间结果
    frameErrorRate    
    
    SNR_ThisPoint = SNR(snr_index);
    
    noiseVar = 10^(-SNR_ThisPoint/10);
    
    hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Es/No)', ...
        'EsNo',SNR_ThisPoint);
    
    hDemod = comm.BPSKDemodulator('DecisionMethod','Log-likelihood ratio', ...
        'Variance',noiseVar);
    
    frameErrorNumber = zeros(1,frameLoop(snr_index));                       %取值为1，表示产生误帧
    
    for frmIdx = 1:frameLoop(snr_index)
        % 生成初始的信息序列，维度 frmLen*1，每个元素随机取值0或者1
        data = randi([0 1],frmLen,1);

        % 进行Turbo码编码
        encodedData = step(hTEnc,data);

        % 进行调制
        modSignal = step(hMod,encodedData);

        % 添加噪声，得到接收信号
        receivedSignal = step(hChan,modSignal);

        %进行解调
        demodSignal = step(hDemod,receivedSignal);

        %进行Turbo译码
        receivedBits = step(hTDec,-demodSignal);

        %进行译码结果判断，统计误比特率
%         errorStats = step(hError,data,receivedBits);
        if(sum(abs(receivedBits-data))~=0)
            frameErrorNumber(frmIdx) = 1;
        end
        
    end
    
    frameErrorRate(snr_index) = mean(frameErrorNumber);    
    
end

%打印结果
frameErrorRate

%绘制FER曲线
semilogy(SNR,frameErrorRate,'*-b')
grid on

xlabel('SNR')
ylabel('FER')




