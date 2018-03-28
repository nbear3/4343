% variable: n double input
% variable: c double input
% variable: R_base double input
% variable: GW double output
% variable: BL double output


%% EC135 Parameters
EW_base = 3208;
GW_base = 6415;
Vtip_base = 41.364*16.75;
A_base = pi*R_base^2;
MCP_base = 1122;
DL = GW_base/A_base;
% Initial Scaling
W_engine_base = engine_weight(MCP_base);
W_drive_base = drive_system_weight(GW_base, MCP_base, DL);

% Constants
phi_struct = (EW_base - W_engine_base - W_drive_base) / GW_base;
ec135.Cd_o = 0.00712; % Cd
ec135.k = 1.15;
ec135.DL = DL;
ec135.sigma = n*c/(pi*R_base);
ec135.Vtip = Vtip_base;

%% Mission Parameters
sfc_base = 0.4;
total_payload = 400 + 1600 + 350;

%% Scaling
GW = 5400;

Rf_avail = 0;
Rf_req = inf;

epsilon = 1e-5;
while (abs(Rf_req - Rf_avail) > epsilon)
    W = GW;

    %% Performance Requirements
    req1 = struct('V', 0, 'V_c', 0, 't', 0, 'h', 13000, 'W', .75*W);
    req2 = struct('V', 0, 'V_c', 0, 't', 0, 'h', 11000, 'W', .90*W);
    req3 = struct('V', 0, 'V_c', 900/60, 't', 0, 'h', 0, 'W', .75*W);
    req4 = struct('V', 165*1.68, 'V_c', 0, 't', 0, 'h', 0, 'W', .8*W);
    req5 = struct('V', 0, 'V_c', 0, 't', 0, 'h', 0, 'W', .8*W);

    [~, ~, P1_sl] = phase_calc(GW, req1.W, 0, ec135, req1);
    [~, ~, P2_sl] = phase_calc(GW, req2.W, 0, ec135, req2);
    [~, ~, P3_sl] = phase_calc(GW, req3.W, 0, ec135, req3);
    [~, ~, P4_sl] = phase_calc(GW, req4.W, 0, ec135, req4);

    P5_req = inf;
    for V=110:150
        req5.V = V;
        [~, ~, P5_req_sl] = phase_calc(GW, req5.W, 0, ec135, req5);
        P5_req = min(P5_req, P5_req_sl);
    end
    P5_sl = GW*(1300/60)/550 + P5_req;

    MCP = max([P1_sl; P2_sl; P3_sl; P4_sl; P5_sl]);

    %% Engine Selection
    W_engine = engine_weight(MCP);
    W_drive = drive_system_weight(GW, MCP, DL);
    Rf_avail = 1-phi_struct-total_payload./GW-(W_engine+W_drive)./GW; 
    gamma = MCP./MCP_base;
    sfc = sfc_base * (-.00932*gamma.^2+.865*gamma+.445)./(gamma+.301);

    %% Mission Phases
    mission(1) = struct('V', 0, 'V_c', 0, 't', 20*60, 'h', 0);
    mission(2) = struct('V', 65*1.68, 'V_c', 350/60, 't', 7500/(350/60), 'h', 7500/2);
    mission(3) = struct('V', 118*1.68, 'V_c', 0, 't', 3600*300/118, 'h', 7500);
    mission(4) = struct('V', 0, 'V_c', 0, 't', 25*60, 'h', 0);
    
    Rf_req = 0;
    for phase = mission 
        F_req = phase_calc(GW, W, sfc, ec135, phase);
        W = W - F_req;
        Rf_req = Rf_req + F_req./GW;
    end
    
    GW = GW + (Rf_req - Rf_avail)*GW;
end

BL = GW/(density(0)*GW/DL*ec135.Vtip^2)/ec135.sigma;

%% Helper Functions
function [F_req, P_req, P_req_sl] = phase_calc(GW, W, sfc, rotor, phase)
    
    Cd_o = rotor.Cd_o;
    k = rotor.k;
    DL = rotor.DL;
    sigma = rotor.sigma;
    Vtip = rotor.Vtip;
    
    A = GW/DL;
    f = -2e-8*GW.^2 + 0.0011*GW + 1.9143; % ft2

    % Phase Parameters
    V = phase.V;
    V_c = phase.V_c;
    hours = phase.t/3600;
    rho = density(phase.h);
    
    % Calculation
    D = .5*rho*V^2*f;
    alpha = atan(D./W);
    T = sqrt(D.^2+W.^2);
    Ct = T./(rho*A*Vtip^2);
    mu = V.*cos(alpha)/Vtip;
    lambda = (V_c/2+sqrt((V_c/2).^2+DL/(2*rho)))/(Vtip);
    Ku = sqrt((-(mu./lambda).^2+sqrt((mu./lambda).^4+4))/2);
    Cp_i = k*Ct.*lambda.*Ku;
    Cp_o = sigma*Cd_o/8*(1+4.65*mu.^2);
    Cp_p = .5*f.*(mu.^3)./A;
    Cp = Cp_i + Cp_o + Cp_p;
    Cp_tail = .05*Cp;
    
    P_req = (Cp+Cp_tail).*(rho*A*(Vtip)^3)/550; % in HP
    P_req_sl = P_req*density(0)/rho;
    
    % Update weight
    F_req = P_req*hours.*sfc;
end

function rho = density(h)
    % Calculate Desnity
    rho = .00238*(1-.00198*h/288.16)^4.2553;
end

function W_engine = engine_weight(MCP) 
    W_engine = (.1054*MCP.^2+358*MCP+2.757e4)./(MCP+1180);
end

function W_drive = drive_system_weight(GW, MCP, DL)
    W_drive =  (525*(GW/1000).^(-1.14))./(((GW./MCP).^.763)*(DL).^.381);
end
