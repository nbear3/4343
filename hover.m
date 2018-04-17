clc
clearvars -except airfoil
close all

c = .8;
theta_t = 5;
tr = .47;
h = 3000*3.28084; % ft

intmesh.OR = 650;
intmesh.t = theta_t;
intmesh.R = 1.3*3.28084;
intmesh.c = c;
intmesh.tr = tr;
intmesh.n = 2;
intmesh.thrust_loss = cosd(13);
intmesh.power_loss = 1.1;
intmesh.download = 0.9;

airfoil.root = xlsread('airfoils/3Rnaca23012.xlsx');
airfoil.tip = xlsread('airfoils/3Tnaca23012.xlsx');
airfoil.ratio = .5;

W = 600*2.2046;
T_req = W/2; % per rotor
V_inf = 0;
max_thr = 25;
epsilon = 1e-3;
[P, P_sl, FM] = thrust_bemt(T_req, epsilon, intmesh, h, V_inf, airfoil, max_thr);

CT_over_CQ = T_req/(P*550)*intmesh.OR; % CT/CQ

% Momentum Theory
Ct = T_req/(density(h)*pi*intmesh.R^2*intmesh.OR^2);
P_momentum = 2*Ct*sqrt(Ct/2)*(density(h)*pi*intmesh.R^2*intmesh.OR^3)/550/FM;

Total_P_req = 2*P_sl;