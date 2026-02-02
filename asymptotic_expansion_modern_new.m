function [R0,R1] = asymptotic_expansion_modern_new(h,v)
kk = 13;
P = zeros(kk,1);
Q = zeros(kk,1);
d = sqrt(h^2+v^2);
alpha = -v/d;
P(1) = 1; P(2) = alpha;
Q(1) = 1; Q(2) = 3*alpha;
for n = 3:kk
    P(n) = (2*(n-1)-1)*alpha*P(n-1)-((n-1)-1)^2*P(n-2);
    Q(n) = (2*(n-1)+1)*alpha*Q(n-1)-((n-1)+1)*((n-1)-1)*Q(n-2);
end
sumP = 0;
sumQ = 0;
for n = 1:kk
    sumP = sumP + P(n)/(d^n);
    sumQ = sumQ + Q(n)/(d^(n+2));
end
R0 = -1*sumP;
R1 = -h*sumQ;

if (v>=-14.5)
    R0 = R0 - pi*exp(v)*bessely(0,h);
    R1 = R1 - pi*exp(v)*bessely(1,h);
end