%��BEC�ŵ��������е�ЧAWGN�ŵ�����,��L��ʾ������ĸ��ŵ�����,I��ʾ˳��%
%����n�Σ��볤N=2^n%
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

        


