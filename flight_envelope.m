clear
clc
close all

S = 5.884*(3.28084)^2;
g = 32.17;
Cl_max = 1.4299;
V = linspace(0, 150/.592484, 500);
Vkts = V*.592484;
[W, intmesh, airfoil] = load_intmesh();

for h = [0 3000*3.28]
    
    rho = density(h);
    L = .5*Cl_max*rho*V.^2*S;
    n_aero = L/W;
    
    j = 1;
    for thr = 34:.2:39
        [T(j), ~, ~, ~] = bemt(intmesh, thr, h, 0, airfoil);
        j=j+1;
    end
    n_rotor = 2*max(T)/W; 
    
    [Vkts_power, n_power] = power_limit(h);    

    %% Plots
    figure
    hold on
    p(1) = plot(Vkts_power, n_power, 'r'); % power limit

    xl = xlim;
    i = InterX([Vkts; n_aero], [[0 100]; [n_rotor n_rotor]]);
    p(2) = plot([xl(1) i(1)], [n_rotor n_rotor], 'g'); % rotor stall
    plot([i(1) xl(2)], [n_rotor n_rotor], 'g--')
    p(3) = plot(Vkts(Vkts >= i(1)), n_aero(Vkts >= i(1)), 'b'); % wing stall
    plot(Vkts(Vkts < i(1)), n_aero(Vkts < i(1)), 'b--')
    i = InterX([Vkts; n_aero], [[0 200]; [2.5 2.5]]);
    p(4) = plot([xl(1) i(1)], [2.5 2.5], 'c'); % structural for rotor
    p(5) = plot(xlim, [4.4 4.4], 'k'); % typical utility category airplane
    plot([i(1) xl(2)], [2.5 2.5], 'c--')
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);

    % title('Flight Envelope')
    xlabel('Speed (kts)')
    ylabel('n')
    legend(p, 'Power',  'Aero (Rotor)', 'Aero (Wing)', 'Structure (Rotor)', 'Structure (Wing)', 'Location', 'NW')
end