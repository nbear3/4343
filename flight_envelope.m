clear
clc
close all

W = 600*2.20462;
rho = density(0*3.28084);
S = 5.884*(3.28084)^2;
g = 32.17;
Cl_max = 1.4299;

V = linspace(0, 250/.592484);
L = .5*Cl_max*rho*V.^2*S;
n = L/W;

plot(V*.592484, n);
hold on

[W, intmesh, airfoil] = load_intmesh();
h = 0;
j = 1;
for thr = 20:.01:40
    [T(j), P(j), ~, ~] = bemt(intmesh, thr, h, 0, airfoil);
    j=j+1;
end

n_power = 2*interp1(P, T, 280/2)/W;
n_aero_limit = 2*max(T)/W;
plot(xlim, [n_power n_power])
plot(xlim, [n_aero_limit n_aero_limit])

% Vstall = sqrt(2*W/(rho*S*Cl_max))*0.592484; % knots
% n_stall = 3.2*(rho*S*Cl_max) * (Vstall/0.592484)^2/(2*W);

% Vmax = 225;
% 
% num = 100; % number of iteration
% V = linspace(0, 250/0.592484, num); % ft/s
% 
% for i = 1:num
%     n(i) = (rho*S*Cl_max)*V(i)^2 / (2*W);
% end


% hold on
% ylim([0,5])
% xlim([0, 250])
plot(xlim, [3 3])
title('Flight Envelope')
xlabel('V_{\infty} (kts)')
ylabel('Load factor, n')
% legend('Structural', 'Location', 'NW')

