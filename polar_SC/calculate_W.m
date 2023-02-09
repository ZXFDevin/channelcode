function  [LLRs]= calculate_W(lamda,phy,n,B,LLRs)
% W0=W0_initialize;
% W1=W1_initialize;
if(lamda==0)
    return 
end
phy_1=floor((phy+1)/2);
%recurse first, if needed%
if (mod(phy,2)==1)
   [LLRs]=calculate_W(lamda-1,phy_1,n,B,LLRs);
end
for beita=0:1:2^(n-lamda)-1
    if (mod(phy,2)==1)
        LLRs(lamda+1,phy+2^lamda*beita)=LLR_f(LLRs(lamda,phy_1+2^(lamda-1)*2*beita),LLRs(lamda,phy_1+2^(lamda-1)*(2*beita+1)));
    else
        bit=B(lamda+1,phy-1+2^lamda*beita);
        LLRs(lamda+1,phy+2^lamda*beita)=LLR_g(LLRs(lamda,phy_1+2^(lamda-1)*2*beita),LLRs(lamda,phy_1+2^(lamda-1)*(2*beita+1)),bit);
        
    end
end
end           
        
    
function LLR1= LLR_f(a,b)
%LLR1=2*atanh(tanh(a/2)*tanh(b/2));
LLR1=sign(a)*sign(b)*min(abs(a),abs(b));
end


function LLR2= LLR_g(a,b,bit)
LLR2=(-1)^(bit)*a+b;
end
