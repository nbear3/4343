clc
clearvars -except airfoil
close all

c = 0.79;
theta_t = 5;
tr = .47;

T_req = 155.98/2; % per rotor
h = 3000*3.28084;

internmesh.OR = 650;
internmesh.t = theta_t;
internmesh.R = 1.3*3.28;
internmesh.c = c;
internmesh.tr = tr;
internmesh.n = 2;
internmesh.thrust_loss = 1;
internmesh.power_loss = 1;
internmesh.download = 0.95;

if ~exist('airfoil', 'var')
    airfoil = [-15.2500000000000,-0.880900000000000,0.134700000000000;-15,-0.902400000000000,0.124710000000000;-14.7500000000000,-0.936100000000000,0.113070000000000;-14.5000000000000,-0.968700000000000,0.102270000000000;-14.2500000000000,-1.00280000000000,0.0913500000000000;-14,-1.03860000000000,0.0799700000000000;-13.7500000000000,-1.07400000000000,0.0676700000000000;-13.5000000000000,-1.10210000000000,0.0552300000000000;-13.2500000000000,-1.12070000000000,0.0470000000000000;-13,-1.13380000000000,0.0419300000000000;-12.7500000000000,-1.14270000000000,0.0384800000000000;-12.5000000000000,-1.14830000000000,0.0358800000000000;-12.2500000000000,-1.14570000000000,0.0335900000000000;-12,-1.13800000000000,0.0315300000000000;-11.7500000000000,-1.12700000000000,0.0296700000000000;-11.5000000000000,-1.11330000000000,0.0279800000000000;-11.2500000000000,-1.09750000000000,0.0264400000000000;-11,-1.07990000000000,0.0250500000000000;-10.7500000000000,-1.06040000000000,0.0238500000000000;-10.5000000000000,-1.04290000000000,0.0222400000000000;-10.2500000000000,-1.02280000000000,0.0209500000000000;-10,-1.00060000000000,0.0199500000000000;-9.75000000000000,-0.977400000000000,0.0190500000000000;-9.50000000000000,-0.953800000000000,0.0181900000000000;-9.25000000000000,-0.929500000000000,0.0173900000000000;-9,-0.905300000000000,0.0165200000000000;-8.75000000000000,-0.880500000000000,0.0157200000000000;-8.50000000000000,-0.855700000000000,0.0149200000000000;-8.25000000000000,-0.829600000000000,0.0144300000000000;-8,-0.803100000000000,0.0140600000000000;-7.75000000000000,-0.776600000000000,0.0136800000000000;-7.50000000000000,-0.749800000000000,0.0133600000000000;-7.25000000000000,-0.723600000000000,0.0129100000000000;-7,-0.696800000000000,0.0126000000000000;-6.75000000000000,-0.670500000000000,0.0122100000000000;-6.50000000000000,-0.643900000000000,0.0119400000000000;-6.25000000000000,-0.617700000000000,0.0116200000000000;-6,-0.591300000000000,0.0115000000000000;-5.75000000000000,-0.566200000000000,0.0112700000000000;-5.50000000000000,-0.541400000000000,0.0111800000000000;-5.25000000000000,-0.517600000000000,0.0109800000000000;-5,-0.493900000000000,0.0108800000000000;-4.75000000000000,-0.471100000000000,0.0106100000000000;-4.50000000000000,-0.448300000000000,0.0104200000000000;-4.25000000000000,-0.425300000000000,0.0103100000000000;-4,-0.392400000000000,0.0100900000000000;-3.75000000000000,-0.354300000000000,0.00979000000000000;-3.50000000000000,-0.314500000000000,0.00961000000000000;-3.25000000000000,-0.271600000000000,0.00926000000000000;-3,-0.230600000000000,0.00898000000000000;-2.75000000000000,-0.196600000000000,0.00872000000000000;-2.50000000000000,-0.163900000000000,0.00841000000000000;-2.25000000000000,-0.136300000000000,0.00801000000000000;-2,-0.111500000000000,0.00756000000000000;-1.75000000000000,-0.0873000000000000,0.00707000000000000;-1.50000000000000,-0.0636000000000000,0.00656000000000000;-1.25000000000000,-0.0401000000000000,0.00618000000000000;-1,-0.0165000000000000,0.00595000000000000;-0.750000000000000,0.00710000000000000,0.00581000000000000;-0.500000000000000,0.0309000000000000,0.00573000000000000;-0.250000000000000,0.0556000000000000,0.00570000000000000;0,0.0795000000000000,0.00569000000000000;0.250000000000000,0.103400000000000,0.00573000000000000;0.500000000000000,0.127000000000000,0.00580000000000000;0.750000000000000,0.150900000000000,0.00589000000000000;1,0.174400000000000,0.00603000000000000;1.25000000000000,0.206000000000000,0.00621000000000000;1.50000000000000,0.246500000000000,0.00654000000000000;1.75000000000000,0.291000000000000,0.00696000000000000;2,0.329300000000000,0.00731000000000000;2.25000000000000,0.353600000000000,0.00753000000000000;2.75000000000000,0.402100000000000,0.00798000000000000;3,0.426500000000000,0.00822000000000000;3.25000000000000,0.451400000000000,0.00847000000000000;3.50000000000000,0.477000000000000,0.00874000000000000;3.75000000000000,0.503000000000000,0.00905000000000000;4,0.529600000000000,0.00935000000000000;4.25000000000000,0.556400000000000,0.00965000000000000;4.50000000000000,0.583500000000000,0.00994000000000000;4.75000000000000,0.610600000000000,0.0102300000000000;5,0.637600000000000,0.0105500000000000;5.25000000000000,0.665000000000000,0.0108100000000000;5.50000000000000,0.692000000000000,0.0111400000000000;5.75000000000000,0.719000000000000,0.0114900000000000;6,0.746300000000000,0.0117500000000000;6.25000000000000,0.773000000000000,0.0121300000000000;6.50000000000000,0.799300000000000,0.0125600000000000;6.75000000000000,0.826400000000000,0.0128500000000000;7,0.853200000000000,0.0131800000000000;7.25000000000000,0.879500000000000,0.0135800000000000;7.50000000000000,0.904700000000000,0.0141800000000000;7.75000000000000,0.932100000000000,0.0143200000000000;8,0.958900000000000,0.0145600000000000;8.25000000000000,0.984900000000000,0.0149200000000000;8.50000000000000,1.01000000000000,0.0154400000000000;8.75000000000000,1.03710000000000,0.0155800000000000;9,1.06380000000000,0.0157700000000000;9.25000000000000,1.08870000000000,0.0162200000000000;9.50000000000000,1.11550000000000,0.0163900000000000;9.75000000000000,1.14220000000000,0.0165400000000000;10,1.16820000000000,0.0168100000000000;10.2500000000000,1.18940000000000,0.0178300000000000;10.5000000000000,1.20810000000000,0.0192200000000000;10.7500000000000,1.22770000000000,0.0204100000000000;11,1.24840000000000,0.0213200000000000;11.2500000000000,1.26770000000000,0.0223700000000000;11.5000000000000,1.28490000000000,0.0235900000000000;11.7500000000000,1.29780000000000,0.0251900000000000;12,1.30920000000000,0.0267700000000000;12.2500000000000,1.32260000000000,0.0280200000000000;12.5000000000000,1.33340000000000,0.0294000000000000;12.7500000000000,1.33900000000000,0.0308800000000000;13,1.34090000000000,0.0325800000000000;13.2500000000000,1.34070000000000,0.0346200000000000;13.5000000000000,1.33820000000000,0.0371200000000000;13.7500000000000,1.33330000000000,0.0402400000000000;14,1.32620000000000,0.0441100000000000;14.2500000000000,1.31650000000000,0.0489600000000000;14.5000000000000,1.30390000000000,0.0550100000000000;14.7500000000000,1.28770000000000,0.0622800000000000;15,1.26790000000000,0.0704000000000000;15.2500000000000,1.24500000000000,0.0789000000000000;15.5000000000000,1.22080000000000,0.0874000000000000];
end

thr = theta_t;
T = 0;
V_forward = 225; % ft/s

while abs(T - T_req) > 1e-4
    [T, P, P_sl, FM] = bemt(internmesh, thr, h, V_forward, airfoil);
    thr = thr + 1/(1+exp((T-T_req)/T_req))-.5;
end

sigma = internmesh.n*internmesh.c/(pi*internmesh.R); 
CT_over_CQ = T/(P*550)*internmesh.OR; % CT/CQ
BL = 1.5*T/(density(h)*pi*internmesh.R^2*internmesh.OR^2)/sigma; % Blade Loading at 0 m

Total_T_req = 2*T_req;
Total_P_req = 2*P_sl;