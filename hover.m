clc
clearvars -except airfoil
close all

c = .6;
theta_t = 10.1;
tr = 1;
h = 0; % ft

intmesh.OR = 700;
intmesh.t = theta_t;
intmesh.R = 1.3*3.28;
intmesh.c = c;
intmesh.tr = tr;
intmesh.n = 2;
intmesh.thrust_loss = cosd(13);
intmesh.power_loss = 1.1;
intmesh.download = 0.9;


airfoil.root = xlsread('sc1095.xlsx');
airfoil.tip = xlsread('sc1095.xlsx');
airfoil.ratio = 1;


W = 600*2.2;
T_req = W/2; % per rotor
V_inf = 0;
max_thr = 25;
epsilon = 1e-3;
[P, P_sl, FM] = thrust_bemt(T_req, epsilon, intmesh, h, V_inf, airfoil, max_thr);

% sigma = intmesh.n*intmesh.c/(pi*intmesh.R); 
% CT_over_CQ = T_req/(P*550)*intmesh.OR; % CT/CQ
% BL = T_req/(density(h)*pi*intmesh.R^2*intmesh.OR^2)/sigma; % Blade Loading at 0 m

% Momentum
Ct = T_req/(density(h)*pi*intmesh.R^2*intmesh.OR^2);
P_momentum = 2*Ct*sqrt(Ct/2)*(density(h)*pi*intmesh.R^2*intmesh.OR^3)/550/FM;

Total_P_req = 2*P_sl;