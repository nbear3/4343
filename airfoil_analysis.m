clc
clear
close all

c = .6;
theta_t = 10.1;
tr = 0.9;

W = 600*2.20462;
h = 0;

intmesh.OR = 700;
intmesh.t = theta_t;
intmesh.R = 1.3*3.28;
intmesh.c = c;
intmesh.tr = tr;
intmesh.n = 2;
intmesh.thrust_loss = cosd(13);
intmesh.power_loss = 1.1;
intmesh.download = 0.9;

airfoils = load_airfoils;

hold on
for i = 1:length(airfoils)
    % Odd is tip; Even is root
    airfoils(i).ratio = 0.6;
    
    j = 1;
    for thr = theta_t:.1:25
        [T(j), P, P_sl, FM] = bemt(intmesh, thr, h, 0, airfoils(i));
        CT_CQ(j) = T(j)/P/550*intmesh.OR;
        j=j+1;
    end
    
    plot(2*T, CT_CQ)
%     P_hover(i) = thrust_bemt(W, 1e-10, intmesh, h, 0, airfoils(i), 25);    
end 

xlabel('Thrust (lbs)')
ylabel('C_T/C_Q')

x = [W W];
plot(x, ylim, '--')
plot(2*x, ylim, '--')
legend('Huey', 'Osprey', 'K-max', 'UH-60A', 'Last Year', 'W', '2W')
title('CT over CQ for Common Airfoil Combinations in Hover at SL')















