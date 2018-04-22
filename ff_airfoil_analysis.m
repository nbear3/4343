clc
clear
close all

h = 3000*3.28084;

ff_data = xlsread('forward_flight.xlsx');
[W, intmesh, airfoil] = load_intmesh();
intmesh.thrust_loss = 1;
intmesh.power_loss = 1;
intmesh.download = .95;

airfoils = load_airfoils();

hold on
V_kts = ff_data(:,2)*0.592484;

for airfoil = airfoils
    airfoil.ratio = 0.5;
    j = 1;
    for r = ff_data'    
        T = r(1)/2;
        V = r(2);
        [~, P, ~] = thrust_bemt(T, 1e-3, intmesh, h, V, airfoil, inf);
        P_req(j) = 2*P;
        j=j+1;
    end
    
    plot(V_kts, P_req)
end

y = [300 300];
plot(xlim, y, '--')

xlabel('V_{\infty} (kts)')
ylabel('Power Required @ SL (HP)')
legend('Naca63-210 (R) & Naca23012 (T)', 'Naca 23012 (R&T)','Naca63-210 (R) & VR-1 (T)', 'Naca23012 (R) & VR-1 (T)', 'RC0864C (R) & Naca23012 (T)')
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
fprintf('\nMax V_{\infty}: %.2f kts\n', V_max);















