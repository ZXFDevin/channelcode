clear all;
clc;
profile on;
% core = 5;
% p = parpool(core);
%============================可调参数=====================================
N = 256;                %码长
K = 128;                %信息比特位数
rate = K/N;             %码率
n = log2(N);
max_iter = 80;
[upnode_M, downnode_M] = get_index(n);

Eb = 1;                    %比特能量归一化
SNR = 2;           %信噪比范围（Eb/N0）
ber = zeros(1, length(SNR));
fer = zeros(1, length(SNR));
ber_c = zeros(1, length(SNR));
fer_c = zeros(1, length(SNR));
iter = zeros(1, length(SNR));
loop_min = 10000;          %最小帧数
loop_max = 1000000;       %最大帧数
fer_err_max = 100;        %最大错误帧数
% 比特翻转规则
index_rule=zeros(1,N);
for index=1:1:N
    b = dec2bin(index-1,n);
    c = dec2bin(index-1,n);
    for i=1:1:length(b)
        c(i) = b(length(b)-i+1);
    end
    index_rule(index) = bin2dec(c)+1;
end


% BEC构造
snr_select = 2;
sigma_1 = sqrt(Eb/(10^(snr_select/10))/2/rate);
HY = quad(@(y) mut_inf(y,sigma_1^2),-20*sigma_1,20*sigma_1)/log(2);
HN = log2(sqrt(2*pi*exp(1))*sigma_1);
C = HY - HN;
P = 1-C;
[L,I] = bec_channel_polarization(P,N);   %得到极化信道可靠性排序
A_C = sort(I(K+1:N),'ascend');  %冻结比特位置信息
A = sort(I(1:K),'ascend');      %信息比特位置信息


G = get_generator(N);

%填充比特
frozen_array = zeros(1,N);
frozen_array(A_C) = 1;
fprintf('开始循环');
%============================主循环=======================================
for loop_SNR = 1:length(SNR)
    snr = 10^(SNR(loop_SNR)/10);            %snr = Eb/N0
    N0 = Eb/snr;
    sigma = sqrt(N0/2/rate);                %N0=2*rate*sigma^2 
    
    ber_err = 0;
    fer_err = 0;
    err = 0;
    tot_iter = 0;
    different_num=0;
    dnf=0;
    ecc_c=0;
     ber_err_c = 0;
    fer_err_c = 0;
    
    % 每一帧循环
    for loop_frame=1:loop_max
        data = randsrc(1,K,[0,1]);     %随机生成K位信息比特
        msg_src = zeros(1,N);
        msg_src(A) = data;             %比特混合
%==============================编码=======================================
        msg_encode = mod(msg_src*G,2);
%============================BPSK调制=====================================
        msg_mapping = (1-2*msg_encode);
%=============================加噪声======================================
        in_AWGN = msg_mapping + sigma*(randn(1,N)+1i*randn(1,N)); 
%==============================解调=======================================
        msg_in = real(in_AWGN);
%==============================译码=======================================
        [msg_decode, tot_iter,flag]=BP_decode(msg_in,N,frozen_array,sigma,max_iter,index_rule,upnode_M, downnode_M,G,tot_iter);

       [code_decided_c,x_e,total_iter,check_pass]=BP_decoder_c(msg_in,sigma,frozen_array,length(frozen_array));
        %========================== ==错误统计====================================
        code_decided=code_decided_c';
        different_num= length(find(code_decided~= msg_decode));
         
         
                 if different_num~=0
            dnf=dnf+1;
        end
        err= length(find(msg_src~= msg_decode));
        if err~=0
          ber_err=ber_err+err;
          fer_err=fer_err+1;
        end
        
        err_c= length(find(msg_src~= code_decided_c'));
        if err_c~=0
          ber_err_c=ber_err_c+err_c;
          fer_err_c=fer_err_c+1;
        end
        
        %打印进度
        if(mod(loop_frame,500)==0)
            fprintf('The SNR:%d\tThe frame number:%d\tThe error number:%d\n',SNR(loop_SNR),loop_frame,fer_err); 
            fprintf('The SNR:%d\tThe frame number:%d\tThe error number:%d\n',SNR(loop_SNR),loop_frame,fer_err_c); 
        end
        %终止条件
        if(fer_err > fer_err_max && loop_frame > loop_min)
            break;
        end
    end
    %统计各信噪比下纠错性能
    ber(loop_SNR) = ber_err/(N*loop_frame);
    fer(loop_SNR) = fer_err/(loop_frame); 
    ber_c(loop_SNR) = ber_err_c/(N*loop_frame);
    fer_c(loop_SNR) = fer_err_c/(loop_frame); 
    iter(loop_SNR) = tot_iter/(loop_frame);
end
figure
semilogy(SNR,ber,'-*b',SNR,fer,'-or');
figure 
semilogy(SNR,ber_c,'-*b',SNR,fer_c,'-or');
profile viewer
