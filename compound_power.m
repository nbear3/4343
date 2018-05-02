function [Vkts, THP, THP_1, THP_2, THP_3] = compound_power(h)
    % Given height h in feet
    %
    % Return:
    %   Vkts: speed in knots
    %   THP: Total Horse Power
    %   THP_1: Parasite Power
    %   THP_2: Shared Transition Power
    %   THP_3: Induced + Profile power
    
    % Gross Weight
    W = 600*2.20462;
    rho_sl = density(0); %slugs/ft^3
    rho = density(h);
    sigmad = rho/rho_sl; %density ratio of air

    %% Compound Helicopter Equations
    e = 0.8; % Oswald's Efficiency Factor
    V = linspace(0, 360); % in ft/s

    S = 5.884*(3.28084)^2;
    Cl_max = 1.45298;
    Vstall = sqrt(2*W/(rho*S*Cl_max)); % Transition point in ft/s
    KL = max(1 - (V./Vstall), 0); % V is speed of helicopter in ft/s
    b = 17.4; %span in ft
    nh = 0.65; %helicopter efficiency, using figure of merit value
    np = 0.65; %propeller efficiency, using figure of merit value

    c = 0.7;
    R = 1.1*3.28; % rotor radius in ft
    A = 2*pi*R^2; % area of both intermeshing rotors in ft^2
    v = sqrt(W/(2*rho*A)); % hover induced velocity

    Ku = sqrt((-1.*(V./v).^2+sqrt((V./v).^4+4))/2); % hover factor
    Vtip = 650; % tip speed in ft/s
    mu = V./Vtip; % normalized speed
    Kmu = 1 + 4.65.*mu.^2; % mu factor

    f = 1.8; % equivalent flat plate area
    Ab = 4*c*R; %blade area
    w = W/A; % disk loading in lbs/ft^2
    B = 0.97; % blade tip loss factor
    Vt0 = Vtip/100; % tip speed in fps/100
    Cd0 = 0.007; %drag coefficient
    IHP = 1.2*W/550*sqrt(w/(2*rho)); % Induced Power
    RHP = sigmad*Ab*Vt0^3*Cd0/1.85; % Profile Power

    THP_1 = ((sigmad.*f.*V.^3)./(146000*1.46667^3*np));
    THP_2 = 0.332./(sigmad*e*np*V).*((1-KL)*W/b).^2;
    THP_2(1) = 0; % account for V=0 at hover
    THP_3 = (KL.^1.1.*Ku.*IHP + Kmu.*RHP);
    THP = THP_1 + THP_2 + THP_3; 

    Vkts = V * 0.592484; % converting ft/s to knots
end
