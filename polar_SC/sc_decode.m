function msg_decoded= sc_decode(in_msg,N,SNR,rate,frozen_array)
%³õÊ¼»¯%
 snr=10^(SNR/10);
 sigma=sqrt(1/(2*rate*snr));             %N0=2*sigma^2  
 n=log2(N);
 LLRs=zeros(n+1,N);
 %W0=zeros(n+1,N);
 %W1=zeros(n+1,N);
 B=zeros(n+1,N);
for i=1:1:N
%     syms y;
%     W0(1,i)=1-int(1/(sqrt(2*pi)*sigma)*exp(-(y+1)^2/(2*sigma^2)),y,in_msg(i),inf);
%     W1(1,i)=1-int(1/(sqrt(2*pi)*sigma)*exp(-(y-1)^2/(2*sigma^2)),y,-inf,in_msg(i));
       %W0(1,i)=1/(sqrt(2*pi)*sigma)*exp(-(in_msg(i)-1)^2/(2*sigma^2));
       %W1(1,i)=1/(sqrt(2*pi)*sigma)*exp(-(in_msg(i)+1)^2/(2*sigma^2));
       LLRs(1,i)=2*in_msg(i)/(sigma^2);
end
%main loop%
for j=1:1:N
    [LLRs]=calculate_W(n,j,n,B,LLRs);
    if (frozen_array(j))
        B(n+1,j)=0;
    else
        if(LLRs(n+1,j)>0)
           B(n+1,j)=0;
        else
           B(n+1,j)=1; 
        end
    end
    if (mod(j,2)==0)
        B=update_B(n,j,n,B);
    end
end
msg_decoded=B(n+1,:);
end
    
    
    
    
