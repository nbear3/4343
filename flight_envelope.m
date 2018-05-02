clear
clc
close all

W = 600*2.20462;
rho = density(0*3.28084);
S = 5.884*(3.28084)^2;
g = 32.17;
Cl_max = 1.4299;

V = linspace(0, 150/.592484, 500);
L = .5*Cl_max*rho*V.^2*S;
n_aero = L/W;

Vkts = V*.592484;
[W, intmesh, airfoil] = load_intmesh();
p = 0;
j = 1;
for thr = 30:.5:36
    [T(j), P(j), ~, ~] = bemt(intmesh, thr, p, 0, airfoil);
    j=j+1;
end

[Vkts_power, n_power] = power_limit();
n_rotor = 2*max(T)/W; 

%% Plots
hold on
p(1) = plot(Vkts_power, n_power, 'r'); % power limit

xl = xlim;
i = InterX([Vkts; n_aero], [[0 100]; [n_rotor n_rotor]]);
p(2) = plot([xl(1) i(1)], [n_rotor n_rotor], 'g'); % rotor stall
plot([i(1) xl(2)], [n_rotor n_rotor], 'g--')
p(3) = plot(Vkts(Vkts >= i(1)), n_aero(Vkts >= i(1)), 'b'); % wing stall
plot(Vkts(Vkts < i(1)), n_aero(Vkts < i(1)), 'b--')

p(4) = plot(xlim, [4.4 4.4], 'k'); % typical utility category airplane
i = InterX([Vkts; n_aero], [[0 200]; [2.5 2.5]]);
p(5) = plot([xl(1) i(1)], [2.5 2.5], 'c'); % structural for rotor
plot([i(1) xl(2)], [2.5 2.5], 'c--')
set(findall(gca, 'Type', 'Line'),'LineWidth',1);

title('Flight Envelope')
xlabel('Speed (kts)')
ylabel('Load factor, n')
legend(p, 'Power',  'Aero (Rotor)', 'Aero (Wing)', 'Structure (Rotor)', 'Structure Wing', 'Location', 'NW')