clc
clear
close all

h = 0*3000*3.28084;

ff_data = xlsread('forward_flight.xlsx');
[W, intmesh, airfoil] = load_intmesh();
intmesh.thrust_loss = 1;
intmesh.power_loss = 1;
intmesh.download = .95;

j = 1;
for r = ff_data'    
    T = r(1)/2;
    V = r(2);
    [~, P, ~] = thrust_bemt(T, 1e-3, intmesh, h, V, airfoil, inf);
    P_req(j) = 2*P;
    j=j+1;
end

V_kts = ff_data(:,2)*0.592484;
plot(V_kts, P_req)
hold on

V_max = interp1(P_req, V_kts, 300);


P_req./V_kts'

xlabel('V_{\infty} (kts)')
ylabel('Power Required @ SL (HP)')

yl = ylim;
y = [300 300];
plot(xlim, y, '--')
plot(V_max, 300, 'o')
plot([V_max V_max], [yl(1) 300], '--')
legend('P Req', 'MCP Avail', 'V Max')
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
fprintf('\nMax V_{\infty}: %.2f kts\n', V_max);















