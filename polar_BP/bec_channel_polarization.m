%以BEC信道擦除进行等效AWGN信道极化,用L表示极化后的各信道容量,I表示顺序%
%极化n次，码长N=2^n%
function [L,I] = bec_channel_polarization(P,N)
n=log2(N);
C=zeros(n,n);
C(1,1)=1-P;
for i=1:1:n
    for j=1:1:2^(i-1)
        C(i+1,2*j-1)=C(i,j)^2;
        C(i+1,2*j)=2*C(i,j)-C(i,j)^2;
    end
end
[L,I] = sort(C(n+1,:),'descend');

        


