
function u = mut_inf(y,sigma2) 

p =(1/(sqrt(2*pi*sigma2)))*exp(-0.5*(y-1).^2/sigma2);
p = p + (1/(sqrt(2*pi*sigma2)))*exp(-0.5*(y+1).^2/sigma2);
p = 0.5*p;
u = -p.*log(p);
return;

