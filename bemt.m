function [T, P, P_sl, FM] = bemt(rotor, thr, h, climb, airfoil)

    sn = 1000; % number of sections

    % Thrust: Total Thrust Produced
    % Power: Total Power Required
    % FM: Figure of Merit

    % rotor.t: Rotor Blade Twist
    % rotor.R: Rotor Blade Radius (float)
    % rotor.c: Rotor Blade Mean Chord (float)
    % rotor.tr: Rotor Blade Taper Ratio (float)
    % rotor.OR: Rotor Blade Tip Speed (Omega*R) (float)
    % rotor.n: Number of blades (float)
    % rotor.power_loss: Power loss multiplier (float)
    % rotor.thrust_loss: Thrust loss multiplier (float)
    % rotor.download: Download Efficiency (0 < download < 1) (float)
    
    % thr: Rotor Blade Root Incidence Angle (float)
    % h: Altitude (either 3000 or 0) (float)
    % climb: Rate of climb ft/s
    
    % airfoil: Airfoil data (struct)
    % airfoil.root: Data for the root airfoil
    % airfoil.tip: Data for the tip airfoil
    % airfoil.ratio: Percent of blade to switch to tip airfoil
    
    t = rotor.t;
    R = rotor.R;
    c = rotor.c;
    tr = rotor.tr;
    OR = rotor.OR; 
    n = rotor.n;
    
    % Density
    rho = density(h);

    %% Retrieve Airfoil Data
    af_cutoff = round(sn*airfoil.ratio);
    
    % Root
    root.aoa = airfoil.root(:, 1); %Angle of Attack
    root.Cl = airfoil.root(:, 2); %Coefficient of Lift
    root.Cd = airfoil.root(:, 3); %Coefficient of Drag
    root.index = 1:af_cutoff-1;
    
    % Tip
    tip.aoa = airfoil.tip(:, 1); 
    tip.Cl = airfoil.tip(:, 2);
    tip.Cd = airfoil.tip(:, 3); 
    tip.index = af_cutoff:sn;
    
    Cla = ones(1, sn);
    for af = [root tip]
        i = af.aoa<5 & af.aoa>-5;
        a = polyfit(af.aoa(i), af.Cl(i), 1); % lift curve slope.
        Cla(af.index) = a(1)*180/pi; %Lift/Curve Slope
    end

    %% Geometry Calculations
    rp = linspace(1/sn, 1, sn); % Radius Proportion
    ru = R * rp; % Radius in Desired Units
    
    cr = 2 * c / (tr + 1); % Chord at Root
    ct = cr * tr; % Chord at Tip
    c = cr - (cr - ct) * rp; % Chord along Span

    tht = thr - t; % Tip Pitch Angle
    th = linspace(thr,tht, sn);

    dA = c*R/sn;
    sigma = n*c./(2*pi*ru); % Elemental Solidity

    % Change of Reference Radius for Center of Element
    rp = rp - .5/sn;
    ru = ru - .5/sn * R;

    omega = OR / R; %Angular Velocity
    Or = omega * ru;
    
    lamc = climb/OR; % Nondimensional Climb Velocity
    lam = sqrt((sigma.*Cla/16-lamc/2).^2+sigma.*Cla/8.*deg2rad(th).*rp)-(sigma.*Cla/16-lamc/2); % Nondimensional Inflow Velocity
    v = lam*OR; % Induced Velocity

    phi = atand(v./Or); % Angle of Incoming Fluid
    alpha = th - phi; % Angle of Attack

    % Coefficient of Lift Lookup
    Clr = ones(1, sn);
    Cdr = ones(1, sn);
    for af = [root tip]
        Clr(af.index) = interp1(af.aoa, af.Cl, alpha(af.index));
        Cdr(af.index) = interp1(af.aoa, af.Cd, alpha(af.index));
        
        % Accounts for Extrapolated Data
        Clr(isnan(Clr)) = 0;
%         [Cl_max, max_i] = max(af.Cl);
%         i = isnan(Clr);
%         Clr(i) = Cl_max-Cla(i).*deg2rad(alpha(i)-af.aoa(max_i));
    end
    Cdr(isnan(Cdr)) = .2; % overestimate

    %% Rotor Performance Calculations
    Vinf = sqrt(Or .^2 + v .^2); % Incoming Flow Velocity
    dL = .5 .* Vinf .^2 .* Clr .* rho .* dA; % Elemental Lift
    dD = .5 .* Vinf .^2 .* Cdr .* rho .* dA; % Elemental Drag

    dT = dL .* cosd(phi) - dD .* sind(phi); % Elemental Thrust
    dPi = ru .* dL .* sind(phi) .* omega; % Elemental Induced Power
    dPo = ru .* dD .* cosd(phi) .* omega; % Elemental Profile Power

    Pi = n*sum(dPi)/550; % Induced Power
    Po = n*sum(dPo)/550; % Profile Power
    P = 1.15*Pi*rotor.power_loss/rotor.download + Po; % Total Power with losses

    T = sum(dT)*n*rotor.thrust_loss; % Total Intermeshing Thrust
    P_sl = P*density(0)/rho; % Power required at sea level
    FM = Pi / P; % Figure of Merit
end