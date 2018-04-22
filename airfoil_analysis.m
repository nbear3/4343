clc
clear
close all

c = .8;
theta_t = 5;
tr = 0.47;

W = 600*2.20462;
h = 1*3000*3.28084;

intmesh.OR = 650;
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
    airfoils(i).ratio = 0.5;
    
    j = 1;
    for thr = theta_t:.1:25
        [T(j), P, P_sl, FM] = bemt(intmesh, thr, h, 0, airfoils(i));
        CT_CQ(j) = T(j)/P/550*intmesh.OR;
        j=j+1;
    end
    T = T/2.20462;
    plot(2*T, CT_CQ)
%     P_hover(i) = thrust_bemt(W, 1e-10, intmesh, h, 0, airfoils(i), 25);    
end 

xlabel('Thrust (kg)')
ylabel('C_T/C_Q')

x = [W/2.20462 W/2.20462];
plot(x, ylim, '--')
plot(1.5*x, ylim, '--')
legend('Naca63-210 (R) & Naca23012 (T)', 'Naca 23012 (R&T)','Naca63-210 (R) & VR-1 (T)', 'Naca23012 (R) & VR-1 (T)', 'RC0864C (R) & Naca23012 (T)', 'W', '1.5W')
%title('CT over CQ for Common Airfoil Combinations in Hover at 3000 m')
title('CT over CQ for Common Airfoil Combinations in Hover at Sea Level')
set(findall(gca, 'Type', 'Line'),'LineWidth',1);















