function [Thrust, Power, Power_sl, FM] = CoaxBEMT(rotor, thr, h, climb)

%Thrust: Total Thrust Produced
%Power: Total Power Required
%FM: Figure of Merit

t = rotor.t;
R = rotor.R;
c = rotor.c;
tr = rotor.tr;
OR = rotor.OR; 
n = rotor.n;
airfoil = rotor.airfoil;

%thr: Rotor Blade Root Incidence Angle (float)
%t: Rotor Blade Twist
%R: Rotor Blade Radius (float)
%c: Rotor Blade Mean Chord (float)
%tr: Rotor Blade Taper Ratio (float)
%OR: Rotor Blade Tip Speed (Omega*R) (float)
%h: Altitude (either 3000 or 0) (float)
%n: Number of Rotors (float)
%climb: Climb Rate (float)
%airfoil: Airfoil Data Excel Sheet (char)

%%Geometry Calculations for Both Rotors
%Retrieve Airfoil Data
num = xlsread(airfoil);
aoa = num(9:end, 1); %Angle of Attack
Cl = num(9:end, 2); %Coefficient of Lift
Cd = num(9:end, 3); %Coefficient of Drag

%Attempts at Lift/Curve Slope Regression
%Cla = (Cl(92) - Cl(73))/5
%Cla = aoa(70:89)\Cl(70:89)
Cla = 2 * pi; %Lift/Curve Slope

rp = .01 : .01 : 1; %Radius Proportion
ru = R * rp; %Radius in Desired Units

cr = 2 * c / (tr + 1); %Chord at Root
ct = cr * tr; %Chord at Tip
c = cr - (cr - ct) * rp; %Chord along Span

tht = thr - t; %Tip Incidence Angle
th = linspace(thr,tht,100);

%Calculation of Elemental Reference Area
dA = zeros(1, 100);
dA(1) = .5 * (cr + c(1)) * .01 * R;
for i = 1 : length(c) - 1
    dA(i+1) = .5 * (c(i) + c(i+1)) *.01 * R;
end

%Calculation of Elemental Disk Area
dDA = zeros(1,100);
dDA(1) = pi * ru(1)^2;
for i = 1 : length(ru) - 1
    dDA(i+1) = pi * (ru(i+1)^2 - ru(i)^2);
end
DA = sum(dDA);

sigma = dA ./ dDA; %Elemental Solidity

%Change of Reference Radius for Center of Element
rp = rp - .005;
ru = ru - .005 * R;

omega = OR / R; %Angular Velocity
Or = omega * ru;

lamc = ones(1, 100) .* climb ./ Or; %Nondimensional Climb Velocity
lam = ((sigma.*Cla./16 - lamc./2).^2 + sigma.*Cla./8. ...
.*(th*pi/180).*rp).^.5    - (sigma.*Cla./16 - lamc./2); %Nondimensional Inflow Velocity
lami = lam - lamc; %Nondimensional Induced Velocity

vi = lami .* OR; %Induced Velocity

phi = atand(vi ./ Or); %Angle of Incoming Fluid
alpha = th - phi; %Angle of Attack

%Density Switch
rho = density(h);

Clr = interp1(aoa, Cl, alpha); %Coefficient of Lift Lookup
Cdr = interp1(aoa, Cd, alpha); %Coefficient of Drag Lookup

%Accounts for Extrapolated Data
for i = 1 : length(Clr)
   if isnan(Clr(i))
       Clr(i) = 0;
   end
   if isnan(Cdr(i))
       Cdr(i) = 0;
   end
end

%%Upper Rotor Performance Calculations
Vinf = (Or .^2 + vi .^2) .^.5; %Incoming Flow Velocity
dL = .5 .* Vinf .^2 .* Clr .* rho .* dA; %Elemental Lift
dD = .5 .* Vinf .^2 .* Cdr .* rho .* dA; %Elemental Drag

dT = dL .* cosd(phi) - dD .* sind(phi); %Elemental Thrust
PipRU = ru .* dL .* sind(phi) .* omega; %Elemental Induced Power
PppRU = ru .* dD .* cosd(phi) .* omega; %Elemental Profile Power
dQ = (PipRU + PppRU) ./ omega; %Elemental Torque

PipRU = sum(PipRU); %Induced Power
PppRU = sum(PppRU); %Profile Power

TpRU = sum(dT); %Thrust Per Rotor
ThrustU = TpRU * n; %Total Thrust
QpRU = sum(dQ); %Torque Per Rotor
PpRU = QpRU * omega / 550; %Power Per Rotor
PowerU = PpRU * n; %Total Power
FMU = PipRU / (PipRU + PppRU); %Figure of Merit

%%Lower Rotor Aerodynamic Calculations
viU = sqrt(ThrustU/2/rho/DA); %Induced Velocity of Upper Rotor
addInf = 2 * viU / OR; %Additional Inflow Velocity of Lower Rotor
reff = sqrt(DA / 2 / pi); %Radius of Lower Rotor Affected by Wake
reffp = ceil(reff / R * 100); %Percentage of Rotor Affected

lamcL = [lamc(1:reffp)+addInf, lamc(reffp+1:end)]; %Lower Rotor Climb Dist.
lam = ((sigma.*Cla./16 - lamcL./2).^2 + sigma.*Cla./8. ...
.*(th*pi/180).*rp).^.5    - (sigma.*Cla./16 - lamcL./2); %Nondimensional Inflow Velocity
lami = lam - lamcL; %Nondimensional Induced Velocity

viL = lami .* OR; %Induced Velocity

phiL = atand(viL ./ Or); %Angle of Incoming Fluid
alphaL = th - phiL; %Angle of Attack

ClrL = interp1(aoa, Cl, alphaL); %Coefficient of Lift Lookup
CdrL = interp1(aoa, Cd, alphaL); %Coefficient of Drag Lookup

for i = 1 : length(Clr)
   if isnan(ClrL(i))
       ClrL(i) = 0;
   end
   if isnan(CdrL(i))
       CdrL(i) = 0;
   end
end

VinfL = (Or .^2 + viL .^2) .^.5; %Incoming Flow Velocity
dLL = .5 .* VinfL .^2 .* ClrL .* rho .* dA; %Elemental Lift
dDL = .5 .* VinfL .^2 .* CdrL .* rho .* dA; %Elemental Drag

dTL = dLL .* cosd(phiL) - dDL .* sind(phiL); %Elemental Thrust
PipRL = ru .* dLL .* sind(phi) .* omega; %Elemental Induced Power
PppRL = ru .* dDL .* cosd(phi) .* omega; %Elemental Profile Power
dQL = (PipRL + PppRL) ./ omega; %Elemental Torque

PipRL = sum(PipRL); %Induced Power
PppRL = sum(PppRL); %Profile Power

TpRL = sum(dTL); %Thrust Per Rotor
ThrustL = TpRL * n; %Total Thrust
QpRL = sum(dQL); %Torque Per Rotor
PpRL = QpRL * omega / 550; %Power Per Rotor
PowerL = PpRL * n; %Total Power
FML = PipRL / (PipRL + PppRL); %Figure of Merit

Thrust = ThrustU + ThrustL; %Total Thrust
Power = PowerU + PowerL; %Total Power
Power_sl = Power*density(0)/rho;
FM = (PipRU + PipRL) / (PipRU + PipRL + PppRU + PppRL); %Figure of Merit
end


function rho = density(h)
    % Calculate Desnity
    rho = .00238*(1-.00198*h/288.16)^4.2553;
end



















