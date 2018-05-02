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

V = V*.592484;


[W, intmesh, airfoil] = load_intmesh();
h = 0;
j = 1;
for thr = 25:.5:36
    [T(j), P(j), ~, ~] = bemt(intmesh, thr, h, 0, airfoil);
    j=j+1;
end

[Vkts, n_power] = power_limit();
n_power_hover = 2*interp1(P, T, 280/2)/W;
n_rotor = 2*max(T)/W;

hold on
plot(V, n_aero, '--');
plot(Vkts, n_power, '--')
plot(xlim, n_power_hover*ones(1, 2), '--')
plot(xlim, [n_rotor n_rotor], '--')
plot(xlim, [3 3], '--')

% i1 = find(n_aero > n_power_hover, 1);
% plot([0 V(i1)], n_power_hover*ones(1, 2), 'k', 'LineWidth', 1.5)
% 
% i2 = find(n_aero > 3, 1);
% plot(V(i1:i2), n_aero(i1:i2), 'k', 'LineWidth', 1.5)

title('Flight Envelope')
xlabel('Speed (kts)')
ylabel('Load factor, n')
legend('Aero (Wing)', 'Power (Compound)', 'Power (Hover)', 'Rotor Stall', 'Structural', 'Location', 'NW')

