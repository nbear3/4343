clear
clc
close all

W = 600*2.20462;
rho = density(0*3.28084);
S = 5.884*(3.28084)^2;
g = 32.17;
Cl_max = 1.4299;
Vstall = sqrt(2*W/(rho*S*Cl_max))*0.592484; % knots
n_stall = 3.2*(rho*S*Cl_max) * (Vstall/0.592484)^2 / (2*W);
Vmax = 225;

num = 100; % number of iteration
V = linspace(0, 250/0.592484, num); % ft/s

for i = 1:num
    n(i) = (rho*S*Cl_max)*V(i)^2 / (2*W);
end


figure (1)
plot(V*0.592484, n)
hold on
plot(xlim, [n_stall, n_stall], '--')
plot([Vmax Vmax], ylim, '--')
title('Flight Envelope')
xlabel('V_{\infty} (kts)')
ylabel('Load factor, n')
legend('Stall Speed Limit', 'Load Factor Limit', 'Dynamic Pressure Limit', 'Location', 'NW')
