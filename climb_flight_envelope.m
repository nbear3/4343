clc
clear
close all

h = 0*3000*3.28084;

[W, intmesh, airfoil] = load_intmesh();

j = 1;
V_climb = 0:5:100;
for V = V_climb    
    [~, P, ~] = thrust_bemt(W/2, 1e-3, intmesh, h, V, airfoil, 25);
    P_req(j) = 2*P;
    j=j+1;
end

plot(V_climb*0.592484, P_req)
hold on

xlabel('V Climb (kts)')
ylabel('Power Required (HP)')

V_max = interp1(P_req, V_climb, 300)*0.592484;


yl = ylim;
y = [300 300];
plot(xlim, y, '--')
plot(V_max, 300, 'o')
plot([V_max V_max], [yl(1) 300], '--')
legend('P Req', 'MCP Avail', 'V Max')
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
fprintf('\nMax V_{\infty}: %.2f kts\n', V_max);













