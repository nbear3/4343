W = 600*2.20462; % lbs
%p = 0.00176420; % slug/ft^3
p = 0.00238
R = 1.3*3.28084; % ft
A = 2*pi*R^2;
n = 4;
Vt = 680;
DL = W/A;
Ct = 1.5*W/(p*A*Vt^2)
sigma = Ct / (0.1);

c = sigma * (2*pi*R) / n
