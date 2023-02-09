% BPSK modulation and LDPC coding
% Parity Check Matrix H is 256x512
% code rate:1/2

clear all;
close all;
clc
tic
K = 256;            %信息位长度
N = 512;            %编码后长度
M = N - K;
R = K/N;            %码率
Eb=2;               %Es=Eb*R=1   
SNRinDB=1:0.5:3;  

ber = zeros(1,length(SNRinDB));                                                 
fer = zeros(1,length(SNRinDB));
loop_min = 1000;    %仿真最小帧数
loop_max = 50000;   %仿真最大帧数
fer_min = 1e-5;     %仿真错误帧数限定
fer_max = 100;

%-------------校验矩阵形式的变化----------------
load 256_512, H1; 
H=sparse(H1);       %将H1中非0元素和其对应的下标存储

hEnc = comm.LDPCEncoder(H);
hDec = comm.LDPCDecoder(H); 

for loop_snr=1:length(SNRinDB)
    loop_snr
    snr=10^(SNRinDB(loop_snr)/10);  %snr in ratio                             
    sigma=sqrt(Eb/(2*snr));         %snr=Eb/N0,N0=2*sigma^2                           
    ber_err = 0;                                                                            
    fer_err = 0;    

    for loop_frame=1:loop_max
        msg_source = logical(randi([0 1], K, 1));  %信源
        
        %----------------LDPC Encode----------------
        msg_coded = step(hEnc, msg_source);   
        
        %----------------BPSK--------------------
        out_map = 1-2*msg_coded;
        
        %-----------------AWGN--------------------
        noise = sigma*(randn(N,1)+i*randn(N,1));
        
        in_msg = real(out_map + noise);   %对于BPSK,虚部均为噪声  
        LLR = 2*in_msg/(sigma^2);         %计算LLR
        
        %--------------LDPC Decode---------------
        msg_decoded = step(hDec,LLR);                      
    
        number0 = sum(abs(msg_decoded(1:K)-msg_source));
        if number0 ~= 0
            ber_err = ber_err + number0;
            fer_err = fer_err + 1;
        end
  
        if (mod(loop_frame,100)==0)
            fprintf('The SNR:%d\tThe frame number:%d\tThe error number:%d\n',SNRinDB(loop_snr), loop_frame,fer_err);   
        end
        if and(fer_err>fer_max, loop_frame>loop_min)
            break;
        end       
    end

    ber(loop_snr) = ber_err/(N*loop_frame);
    fer(loop_snr) = fer_err/(loop_frame);    

    if fer(loop_snr) <fer_min
        break;
    end
    
end 
toc
% % plot results
semilogy(SNRinDB,ber,'b',SNRinDB,fer,'r');
grid on       %这里SNR定义为BPSK符号中每个bit的Eb/N0
