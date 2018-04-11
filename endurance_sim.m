% variable: n double input
% variable: c double input
% variable: R_base double input
% variable: GW double output
% variable: BL double output

clc
clear

% rotor parameter
intmesh.OR = 725;
intmesh.t = theta_t;
intmesh.R = 1.3*3.28;
intmesh.c = c;
intmesh.tr = tr;
intmesh.n = 2;
intmesh.thrust_loss = cosd(13);
intmesh.power_loss = 1.1;
intmesh.download = 0.9;


% Constants 
GW_base = 600*2.2046;
A_base = pi*R_base^2;

% coax.DL = GW_base/(A_base*2);
% % coax.sigma = 2n*c/(pi*R_base);
% coax.OR = 685;
% % coax.losses = 0.1;
% coax.t = -9;
% coax.R = R_base;
% coax.c = c;
% coax.tr = .8;
% coax.n = n;
% coax.airfoil = 'sc1095.xlsx';

[t, MCP] = hover_time(coax);
fprintf('Coax\n')
fprintf('Hover Time 50%%:   %.2f s\n', t);
fprintf('P Required:        %.2f HP\n', MCP);

function [t, MCP] = hover_time(rotor)

    % Parameters
    EW = 400*2.2046;
    GW = 600*2.2046;
    sfc_base = 0.4;
    MCP_base = 430;

    % Initial Scaling
    W_engine_base = engine_weight(MCP_base);
    W_drive_base = drive_system_weight(GW, MCP_base, rotor.DL);
    phi_struct = (EW - W_engine_base - W_drive_base) / GW;

    % Mission Parameters
    total_payload = 100*2.2046;

    % Performance Requirements
    time_step = 1000;
    hover = struct('V', 0, 'V_c', 0, 't', time_step, 'h', 3000*3.28, 'W', GW);
    [~, ~, P1_sl] = phase_calc(GW, hover.W, 0, rotor, hover);
    MCP = P1_sl;

    % Engine Selection
    W_engine = engine_weight(MCP);
    W_drive = drive_system_weight(GW, MCP, rotor.DL);
    Rf_avail = 1-phi_struct-total_payload./GW-(W_engine+W_drive)./GW; 
    gamma = MCP./MCP_base;
    sfc = sfc_base * (-.00932*gamma.^2+.865*gamma+.445)./(gamma+.301);

%     BL = GW/(density(0)*GW/rotor.DL*rotor.Vtip^2)/rotor.sigma;
    
    t = 0;
    W = GW;    
    Rf_50 = Rf_avail/2;
    while Rf_avail > Rf_50
        F_req = phase_calc(GW, W, sfc, rotor, hover);
        W = W - F_req;
        Rf_avail = Rf_avail - F_req/GW;
        t = t + time_step;
    end
end

%% Helper Functions
function [F_req, P_req, P_req_sl] = phase_calc(GW, W, sfc, rotor, phase)
    
%     Cd_o = 0.008;
%     k = 1.15;
%     DL = rotor.DL;
%     
%     A = GW/DL;
%     f = -2e-8*GW.^2 + 0.0011*GW + 1.9143; % ft2

    % Phase Parameters
    V = phase.V;
    hours = phase.t/3600;
    
    Thrust = 0;
    thr = 15;
    while Thrust < W
        [Thrust, P_req, P_req_sl, FM] = CoaxBEMT(rotor, thr, phase.h, phase.V_c);
        thr = thr + .1;
       
    end
    
    % Update weight
    F_req = P_req*hours.*sfc;
end

function W_engine = engine_weight(MCP) 
    W_engine = (.1054*MCP.^2+358*MCP+2.757e4)./(MCP+1180);
end

function W_drive = drive_system_weight(GW, MCP, DL)
    W_drive =  (525*(GW/1000).^(-1.14))./(((GW./MCP).^.763)*(DL).^.381);
end
