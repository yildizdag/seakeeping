function [G,Gx] = freeSurfGreenFuncModern(x1,x2,x3,xi1,xi2,xi3,f)
% X: Field
% Xi: Source

%Parameters:
rho = sqrt((x1-xi1)^2+(x2-xi2)^2);
r = sqrt(rho^2+(x3-xi3)^2);
rp = sqrt(rho^2+(x3+xi3)^2);
h = f*rho;
v = f*(x3+xi3);
d = f*rp;
%Loop over 5 Different Regions:
if (h^2+0.64*v^2) >= 144
    [R0,R1]=asymptotic_expansion_modern_new(h,v);
elseif ((v>=-2.0) && (h<=1.2-0.1*v))
    [R0,R1]=ascending_series(h,v);
elseif h<=-0.7*v
    [R0,R1]=taylor_series_haxis_modern(h,v);
elseif v>=-0.4*h
    [R0,R1] = taylor_series_vaxis(h,v);
else
    [R0,R1] = haskind_representation(h,v);
end
G = (-1/(4*pi))*(1/r+1/rp)-(2*f/(4*pi))*(R0+1i*pi*besselj(0,h)*exp(v));
Grho = (1/(4*pi))*(rho/(r^3)+rho/(rp^3))+(2*f^2/(4*pi))*(R1+1i*pi*besselj(1,h)*exp(v));
Gx1 = Grho*(x1-xi1)/rho;
Gx2 = Grho*(x2-xi2)/rho;
Gx3 = (1/(4*pi))*((x3-xi3)/r^3+(x3+xi3)/rp^3)-(2*f^2/(4*pi))*((1/d)+R0+1i*pi*besselj(0,h)*exp(v));
Gx = [Gx1; Gx2; Gx3];