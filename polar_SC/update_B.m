function  B=update_B(lamda,phy,n,B)
phy_1=floor((phy+1)/2);
for beita=0:1:2^(n-lamda)-1
    B(lamda,phy_1+2^(lamda-1)*2*beita)=mod(B(lamda+1,phy-1+2^lamda*beita)+B(lamda+1,phy+2^lamda*beita),2);
    B(lamda,phy_1+2^(lamda-1)*(2*beita+1))=B(lamda+1,phy+2^lamda*beita);
end
if (mod(phy_1,2)==0)
   B=update_B(lamda-1,phy_1,n,B);
end

    
    
    



