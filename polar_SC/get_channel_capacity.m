% trasnlate the AWGN channel to BEC channel
% clear all;
function [C,P]= get_channel_capacity(SNR)
% SNR=-20:1:20;
% C=zeros(1,length);
snr=10^(SNR/10);
% C(i)=1-2^(-snr);  %近似计算
%精确计算BI-AWGN信道容量%
syms k;
C=1+log2(exp(1))*(-((2*snr/pi)^(1/2))*exp(-snr/2)+(2*snr-1)*qfunc(snr^(1/2))+vpa(symsum(((-1)^k/(k*(k+1))*qfunc(real((snr^(1/2))*(2*k+1)))*exp(2*snr*k*(k+1))),1,Inf)));
C=double(C);
% end
% plot(SNR,C);

%找出对应的BEC信道erasure probability p%
% p=zeros(1,length);
% for i=1:length
    P=1-C;
% plot(SNR,p);

