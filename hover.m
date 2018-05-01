clc
clear
close all

h = 0*3000*3.28084; % ft
[W, intmesh, airfoil] = load_intmesh();
max_thr = 30;
epsilon = 1e-3;
[~, P_sl, ~, dT, ru] = thrust_bemt(W/2, epsilon, intmesh, h, 0, airfoil, max_thr);

% Best Range Speed
T = 129.53;
Vf = 93.67*1.68781;
ff_intmesh = intmesh;
ff_intmesh.thrust_loss = 1;
ff_intmesh.power_loss = 1;
ff_intmesh.download = .95;
[~, P_ff_sl, ~, dT_ff, ru_ff] = thrust_bemt(T/2, epsilon, ff_intmesh, h, Vf, airfoil, 50);

% Momentum Theory
MCP_req = 2*P_sl;

plot(ru/intmesh.R, dT*4.44822)
hold on
plot(ru/intmesh.R, dT_ff*4.44822)

xlabel('r/R')
ylabel('dT (N)')