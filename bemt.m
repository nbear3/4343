function [T, P, P_sl, FM] = bemt(rotor, thr, h, climb, airfoil)
    %Thrust: Total Thrust Produced
    %Power: Total Power Required
    %FM: Figure of Merit

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
    % airfoil: Airfoil Data Excel Sheet (char)
    
    t = rotor.t;
    R = rotor.R;
    c = rotor.c;
    tr = rotor.tr;
    OR = rotor.OR; 
    n = rotor.n;
    
    % Density
    rho = density(h);

    %%Geometry Calculations for Both Rotors
    % Retrieve Airfoil Data
    aoa = airfoil(2:end, 1); %Angle of Attack
    Cl = airfoil(2:end, 2); %Coefficient of Lift
    Cd = airfoil(2:end, 3); %Coefficient of Drag

    % Lift/Curve Slope Regression
    a = polyfit(aoa(aoa<10 & aoa>-10), Cl(aoa<10 & aoa>-10), 1); % lift curve slope.
    Cla = a(1)*180/pi; %Lift/Curve Slope

    rp = .01 : .01 : 1; % Radius Proportion
    ru = R * rp; % Radius in Desired Units

    cr = 2 * c / (tr + 1); % Chord at Root
    ct = cr * tr; % Chord at Tip
    c = cr - (cr - ct) * rp; % Chord along Span

    tht = thr - t; % Tip Incidence Angle
    th = linspace(thr,tht,100);

    % Calculation of Elemental Reference Area
    dA = zeros(1, 100);
    dA(1) = .5 * (cr + c(1)) * .01 * R;
    for i = 1 : length(c) - 1
        dA(i+1) = .5 * (c(i) + c(i+1)) *.01 * R;
    end

    % Calculation of Elemental Disk Area
    dDA = zeros(1,100);
    dDA(1) = pi * ru(1)^2;
    for i = 1 : length(ru) - 1
        dDA(i+1) = pi * (ru(i+1)^2 - ru(i)^2);
    end
    DA = sum(dDA);

    sigma = dA ./ dDA; % Elemental Solidity

    % Change of Reference Radius for Center of Element
    rp = rp - .005;
    ru = ru - .005 * R;

    omega = OR / R; %Angular Velocity
    Or = omega * ru;
    
    lamc = ones(1, 100) .* climb ./ OR; % Nondimensional Climb Velocity
    lam = ((sigma.*Cla./16 - lamc./2).^2 + sigma.*Cla./8. ...
    .*(th*pi/180).*rp).^.5 - (sigma.*Cla./16 - lamc./2); % Nondimensional Inflow Velocity
    lami = lam - lamc; % Nondimensional Induced Velocity
    vi = lami .* OR; % Induced Velocity

    phi = atand(vi ./ Or); % Angle of Incoming Fluid
    alpha = th - phi; % Angle of Attack

    Clr = interp1(aoa, Cl, alpha); % Coefficient of Lift Lookup
    Cdr = interp1(aoa, Cd, alpha); % Coefficient of Drag Lookup

    % Accounts for Extrapolated Data
    for i = 1 : length(Clr)
       if isnan(Clr(i))
           Clr(i) = 0;
       end
       if isnan(Cdr(i))
           Cdr(i) = 0;
       end
    end

    %% Rotor Performance Calculations
    Vinf = sqrt(Or .^2 + vi .^2); % Incoming Flow Velocity
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









