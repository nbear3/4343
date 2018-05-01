clc
clear
close all

[W, intmesh, airfoil] = load_intmesh();
hold on

h = 0;
j = 1;
for thr = 20:.1:25
    [T(j), P, P_sl(j), FM] = bemt(intmesh, thr, h, 0, airfoil);
    j=j+1;
end

T = T/2.20462; % convert to kg
plot(2*T, P_sl*2)

xlabel('Thrust (kg)')
ylabel('Power (HP)')

W = W/2.20462; % weight in kg
P_hover = interp1(2*T, 2*P_sl, W);
xl = xlim;
yl = ylim;
h1 = plot(W, P_hover, 'o');
plot([W W], [yl(1) P_hover], '--')
plot([xl(1), W], [P_hover P_hover], '--')

legend(h1, 'P Req in Hover');
set(findall(gca, 'Type', 'Line'),'LineWidth',1);















