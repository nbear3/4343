% variable: n double input
% variable: c double input
% variable: R_base double input
% variable: GW double output
% variable: BL double output

clc
clear

n = 2; % per rotor
c = 0.5;
R_base = 4.5; 

% Constants 
GW_base = 600*2.2046;
A_base = pi*R_base^2;

coax.DL = GW_base/(A_base*2);
coax.sigma = 2*n*c/(pi*R_base);
coax.Vtip = 685;
coax.losses = 0.1;

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
    time_step = .1;
    hover = struct('V', 0, 'V_c', 0, 't', time_step, 'h', 3000*3.28, 'W', GW);
    [~, ~, P1_sl] = phase_calc(GW, hover.W, 0, rotor, hover);
    MCP = P1_sl;

    % Engine Selection
    W_engine = engine_weight(MCP);
    W_drive = drive_system_weight(GW, MCP, rotor.DL);
    Rf_avail = 1-phi_struct-total_payload./GW-(W_engine+W_drive)./GW; 
    gamma = MCP./MCP_base;
    sfc = sfc_base * (-.00932*gamma.^2+.865*gamma+.445)./(gamma+.301);

    BL = GW/(density(0)*GW/rotor.DL*rotor.Vtip^2)/rotor.sigma;
    
    t = 0;
    W = GW;    
    Rf_50 = Rf_avail/2;
    while Rf_avail > Rf_50
        hover.W = W;
        F_req = phase_calc(GW, W, sfc, rotor, hover);
        W = W - F_req;
        Rf_avail = Rf_avail - F_req/GW;
        t = t + time_step;
    end
end

%% Helper Functions
function [F_req, P_req, P_req_sl] = phase_calc(GW, W, sfc, rotor, phase)
    
    Cd_o = 0.008;
    k = 1.15;
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
    Cp_tail = Cp * rotor.losses;
    
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
