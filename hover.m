clc
clearvars -except airfoil
close all


h = 0*3000*3.28084; % ft

[W, intmesh, airfoil] = load_intmesh();
intmesh.thrust_loss = 1;
intmesh.power_loss = 1;
intmesh.download = .95;

T_req = W/2; % per rotor
V_inf = 0;
max_thr = 25;
epsilon = 1e-3;
[P, P_sl, FM] = thrust_bemt(T_req, epsilon, intmesh, h, V_inf, airfoil, max_thr);

CT_over_CQ = T_req/(P*550)*intmesh.OR; % CT/CQ

% Momentum Theory
Ct = T_req/(density(h)*pi*intmesh.R^2*intmesh.OR^2);
P_momentum = 2*Ct*sqrt(Ct/2)*(density(h)*pi*intmesh.R^2*intmesh.OR^3)/550/FM;

MCP_req = 2*P_sl;