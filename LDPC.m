% BPSK modulation and LDPC coding
% Parity Check Matrix H is 256x512
% code rate:1/2

clear all;
close all;
clc
tic
K = 256;            %��Ϣλ����
N = 512;            %����󳤶�
M = N - K;
R = K/N;            %����
Eb=2;               %Es=Eb*R=1   
SNRinDB=1:0.5:3;  

ber = zeros(1,length(SNRinDB));                                                 
fer = zeros(1,length(SNRinDB));
loop_min = 1000;    %������С֡��
loop_max = 50000;   %�������֡��
fer_min = 1e-5;     %�������֡���޶�
fer_max = 100;

%-------------У�������ʽ�ı仯----------------
load 256_512, H1; 
H=sparse(H1);       %��H1�з�0Ԫ�غ����Ӧ���±�洢

hEnc = comm.LDPCEncoder(H);
hDec = comm.LDPCDecoder(H); 

for loop_snr=1:length(SNRinDB)
    loop_snr
    snr=10^(SNRinDB(loop_snr)/10);  %snr in ratio                             
    sigma=sqrt(Eb/(2*snr));         %snr=Eb/N0,N0=2*sigma^2                           
    ber_err = 0;                                                                            
    fer_err = 0;    

    for loop_frame=1:loop_max
        msg_source = logical(randi([0 1], K, 1));  %��Դ
        
        %----------------LDPC Encode----------------
        msg_coded = step(hEnc, msg_source);   
        
        %----------------BPSK--------------------
        out_map = 1-2*msg_coded;
        
        %-----------------AWGN--------------------
        noise = sigma*(randn(N,1)+i*randn(N,1));
        
        in_msg = real(out_map + noise);   %����BPSK,�鲿��Ϊ����  
        LLR = 2*in_msg/(sigma^2);         %����LLR
        
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
grid on       %����SNR����ΪBPSK������ÿ��bit��Eb/N0
