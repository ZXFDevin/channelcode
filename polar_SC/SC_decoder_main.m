%极化码仿真主程序，包括encode以及decode,在AWGN信道下%
%Communication Letters Systematic Polar Coding%
clear all;
clc;
%=============================可以调整的参数================================

K=128;                 %编码前的数据量，
N=256;                 %总共的子载波数量,码长
n=log2(N);           %block size%
M=N-K;
rate = K/N;           %码率

Ebpsk=1; %energy per bit
Eb=1;
SNRinDB=0:0.5:2;            %这里SNR定义为QPSK符号中每个bit的Eb/N0
ber = zeros(1,length(SNRinDB));
fer = zeros(1,length(SNRinDB));
loop_loop_min = 1000;             %仿真最小帧数
loop_loop_max=10000000;           %仿真最大帧数
fer_err_max = 100;                %仿真最大错误比特数
fer_min = 1e-5;    

G_generator= get_generator(N);
for loop_snr=1:length(SNRinDB)
    snr=10^(SNRinDB(loop_snr)/10);           %snr in ratio.snr=Eb/N0;
    sigma=sqrt(Eb/(2*rate*snr));             %N0=2*sigma^2  
    %对每一个信噪比进行信道极化，并且选定frozen bit的位置%
    HY = quad(@(y) mut_inf(y,sigma^2),-20*sigma,20*sigma)/log(2);
    HN = log2(sqrt(2*pi*exp(1))*sigma);
    C = HY - HN;
    P = 1-C;
    [L,I] = bec_channel_polarization(P,N);
     %用frozen_position指代frozen bits的位置%
    frozen_positions=sort(I(K+1:N),'ascend');
    free_positions = sort(I(1:K),'ascend');
    frozen_positions = sort(frozen_positions,'ascend');
    frozen_array=zeros(1,N);
    frozen_array(frozen_positions)=1;
    err=0;                              %统计译码错误数据   
    ber_err = 0;
    fer_err = 0;
    for loop_frame=1:loop_loop_max    
    data_in=randsrc(K,1,[0 1]); % random data sequence
    msg_source = zeros(N,1);
    msg_source(free_positions) = data_in; % set frozen positions to zero
    %----------------极化码编码----------------%
    msg_coded= encode(data_in,G_generator,frozen_positions,free_positions);
%     msg_coded=rvsl(msg_coded);
    %---------------BPSK mapping-----------------%    
    msg_mapping=(1-2*msg_coded);          
    
%     %---------------AWGN信道-------------------%
     N0=sigma*(randn(N,1)+1i*randn(N,1));  %complex AWGN    
     in_AWGN=msg_mapping+N0;                                   %add AWGN noise
 
     
    %----------------demapping---------------%
    in_msg=real(in_AWGN);  %in BPSK , the receiver just process real part, imag part is all noise.
    
    %-----------------极化码解码-----------------%
     msg_decoded= sc_decode(in_msg,N,SNRinDB(loop_snr),rate,frozen_array);
     err= length(find(msg_decoded~= msg_source'));


       if err~=0
          ber_err=ber_err+err;
          fer_err=fer_err+1;
       end
    
       if (mod(loop_frame,100)==0)
          fprintf('The SNR:%d\tThe frame number:%d\tThe error number:%d\n',SNRinDB(loop_snr),loop_frame,fer_err);   
       end
       if and(fer_err >fer_err_max, loop_frame>loop_loop_min)
          break;
       end
    end
    
    ber(loop_snr) = ber_err/(N*loop_frame);
    fer(loop_snr) = fer_err/(loop_frame); 
    fid_err=fopen('ber_fer.txt','a');
    fprintf(fid_err,'SNR:%6.8f\tBER:%6.8f\tFER:%6.8f\n',SNRinDB(loop_snr),ber(loop_snr),fer(loop_snr));
    fclose(fid_err);
    if fer(loop_snr) <fer_min
        break;
    end
end
% plot results
semilogy(SNRinDB,ber,'b',SNRinDB,fer,'r');
grid on       %这里SNR定义为BPSK符号中每个bit的Eb/N0
