clc
clear 

[GW, rotor, airfoil] = load_intmesh();

phi = 0.55; % includes payload
sfc = 0.578;

% Mission Parameters
payload = 100*2.2046; % kg

% Performance Requirements
h = 1*3000*3.28084; % height (ft)
V = 0; % forward speed (ft/s)

% Fuel Fraction Available
Rf_avail = 1-phi-payload./GW; 

W = GW;    
Rf_50 = Rf_avail/2;
Rf_next = 0;
epsilon = 1e-3;
max_thr = 30; % maximum collective pitch

t = 0;
time_step = 500; % in seconds

i = 1;
while (Rf_avail - Rf_next) > Rf_50
    [P_req, ~, ~] = thrust_bemt(W/2, epsilon, rotor, h, V, airfoil, max_thr);
    F_req = 2*P_req*time_step*sfc/3600; % requires 2*P_req for both rotor
    W = W - F_req;

    Rf_next = F_req/GW;
    Rf_avail = Rf_avail - Rf_next;
    t = t + time_step;  
    
    F_ref(i) = F_req;
    P_ref(i) = 2*P_req;
    t_ref(i) = t;
    i = i + 1;
end

F_left = (Rf_avail - Rf_50)*GW;
[P_req, ~, ~] = thrust_bemt(W/2, epsilon, rotor, h, V, airfoil, max_thr);
t = t + F_left/(2*P_req)/sfc*3600;
h = floor(t/3600);
m = floor(mod(t, 3600)/60);
fprintf('Intermesh Hover Time @ 50%% Fuel:   %d hr %d min\n', h, m);