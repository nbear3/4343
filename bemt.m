% CBEMT with Tip Loss for OH-58D
% Sylvester Ashok

clc

clearvars -except airfoil

if ~exist('airfoil', 'var')
    airfoil = xlsread('sc1095');
end

i=1;
c = 1.1:.1:5;
for c_i = c
    [PL(i), FM(i)] = bemta(c_i, deg2rad(-10), airfoil);
    i = i+1;
end
plot(c, PL)

function [PL, FM] = bemta(c, theta_tw, airfoil)

% OH-58D data
W = 600*2.20462; % lbs
R = 1.3*3.28; % ft
A = 2*pi*R^2; % ft2

root_cutoff = 0.1;
VT = 680;
omega = VT/R;

% rpm = 395;
% % omega = 450*2*pi/60; % rad/sec
% VT = omega*R; % Tip speed ft/sec


aoa = deg2rad(airfoil(2:end, 1)); %Angle of Attack
Cl = airfoil(2:end, 2); %Coefficient of Lift
Cd = airfoil(2:end, 3); %Coefficient of Drag

% c = 1.1; % mean chord ft
Nb = 4; %4 blades
sigma = Nb*c*R/A; % solidity

rho = 0.00176; % slugs/ft3
theta_0 = deg2rad(10); % initial collective
% theta_tw = deg2rad(-5); %2 per foot
taper = 0.7; %tip is 70% of root

a = polyfit(aoa(aoa<10 & aoa>-10), Cl(aoa<10 & aoa>-10), 1); % lift curve slope.
a = a(1);

CLmax = max(Cl);
CD0 = 0.08; % avg drag coeff

% blade sections
n = 100; % 100 blade elements
r = linspace(R/n,R,n+1); % r value
dr = R/n;
m = n*root_cutoff; % number of root cut off sections to ignore
c = linspace(c*2/(1+taper),taper*c,n+1); % sectional chord

% additional loop for trimming collective pitch
o = 30; % number of increments  
T = zeros(o, 1);
P = zeros(o, 1);
Q = zeros(o, 1);
theta_0(1) = theta_0;

for k = 1:o
    
    % initializing vectors
    lambda = zeros(n,1);
    F = ones(n,1);
    v = zeros(n,1);
    U = zeros(n,1);
    dL = zeros(n,1);
    dD = zeros(n,1);
    alpha = zeros(n,1);
    theta = zeros(n,1);
    phi = zeros(n,1);
    dQ = zeros(n,1);
    dT = zeros(n,1);
    
    % Loop 1 for Blade Elements
    for i = m:n
        j = 1;
        l(1) = 0;
        theta(i) = theta_0(k) + theta_tw*(r(i)/R);
        
        %Loop 2 for Tip Loss
        while (true)
            j = j+1;
            l(j) = sigma*a/16/F(i)*(sqrt(1+32*F(i)/sigma/a*theta(i)*r(i)/R)-1);
            f = Nb/2*(1-r(i)/R)/l(j);
            F(i) = 2/pi*acos(exp(-f));
            if abs(l(j)-l(j-1))/l(j-1) <= 0.01
                lambda(i) = l(j);
                l = 0;
                break
            end
        end
        
        % Blade Element Calculations
        v(i) = lambda(i)*VT;
        phi(i) = atan(v(i)/(omega*r(i)));
        alpha(i) = theta(i) - phi(i);
        % stall check
        
        CL = interp1(aoa, Cl, alpha(i)); %Coefficient of Lift Lookup;
        if isnan(CL)
            CL = CLmax - a*(alpha(i)-deg2rad(12));
            CD = .008;
           
        else
            CD = .008; %Coefficient of Lift Lookup;
        end
       
        
        U(i)= sqrt((omega*r(i))^2 + v(i)^2);
        dL(i) = 0.5*rho*U(i)^2*dr*c(i)*CL;
        dD(i) = 0.5*rho*U(i)^2*dr*c(i)*CD;
        dT(i) = Nb*(dL(i)*cos(phi(i)) - dD(i)*sin(phi(i)));
        dQ(i) = Nb*r(i)*(dL(i)*sin(phi(i)) + dD(i)*cos(phi(i)));
        
        
    end
    
    
    T(k) = sum(dT);
    Q(k) = sum(dQ);
    P(k) = Q(k)*omega/550;
    
    theta_0(k+1) = theta_0(k) + deg2rad(0.5);
end

intermesh_loss = 1.1;
CP = P.*550./(rho*A*VT^3)*intermesh_loss;
CT = T./(rho*A*VT^2);
FM = CT.*sqrt(CT./2)./CP;

CT_o = 1.5*W/(rho*A*VT^2);
FM = interp1(CT, FM, CT_o); %Coefficient of Lift Lookup;
CP_o = interp1(CT, CP, CT_o); %Coefficient of Lift Lookup;
PL = CT_o/CP_o;

end

% figure(2)
% plot(r(m:n)./R, dT(m:n))
% xlabel ('nondimensional radius')
% ylabel ('dT')
% figure(3)
% plot(r(m:n)./R, dQ(m:n))
% xlabel ('nondimensional radius')
% ylabel ('dQ')
% figure(4)
% plot(r(m:n)./R, lambda(m:n))
% xlabel ('nondimensional radius')
% ylabel ('lambda')
% figure (1)
% plot(CP,CT)
% xlabel ('CP')
% ylabel ('CT')
% figure (2)
% plot(CT,FM)
% xlabel ('CT')
% ylabel ('FM')








