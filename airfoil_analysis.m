clc
clear
close all

c = 1;
theta_t = 10.1;
tr = 0.9;

W = 600*2.2/2;
h = 3000*3.28;

intmesh.OR = 700;
intmesh.t = theta_t;
intmesh.R = 1.3*3.28;
intmesh.c = c;
intmesh.tr = tr;
intmesh.n = 2;
intmesh.thrust_loss = cosd(13);
intmesh.power_loss = 1.1;
intmesh.download = 0.9;

airfoils = load_airfoils;

for i = 1:length(airfoils)
    % Odd is tip; Even is root
    airfoils(i).ratio = 0.5;
    P_hover(i) = thrust_bemt(W, 1e-10, intmesh, h, 0, airfoils(i), 25);    
end




















