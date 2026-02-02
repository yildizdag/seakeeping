function [R0,R1] = taylor_series_haxis_modern(h,v)
P = zeros(15,1);
P(1) = 1-v;
for n = 2:15
    P(n) = factorial(2*n-2)*(2*n-1-v)/factorial(n)^2+(v^2/n^2)*P(n-1);
end
sum0 = 0;
for n = 1:15
    sum0 = sum0 + P(n)*(-h^2/(4*v^2))^n;
end
sum1 = 0;
for n = 1:15
    sum1 = sum1 + n*P(n)*(-h^2/(4*v^2))^(n-1);
end
R0 = besselj(0,h)*(exp(v)*real(expint(v+1i*0)))+sum0;
R1 = besselj(1,h)*(exp(v)*real(expint(v+1i*0)))+(h/(2*v^2))*sum1;