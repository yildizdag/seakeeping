function phiI = incident(r,freq,beta,zeta)
k = freq^2/9.81; 
phiI = 1i*(9.81*zeta/freq)*exp(k*r(3))*exp(1i*k*(r(1)*cosd(beta)+r(2)*cosd(beta)));
