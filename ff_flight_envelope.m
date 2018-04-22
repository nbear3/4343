clc
clear
close all

%h = 3000*3.28084;
h = [0, 3000*3.28084];
fn.name = {'forward_flight.xlsx', 'forward_flight_alt.xlsx'};

figure (1)

for n = 1:length(h)
    %ff_data = xlsread('forward_flight.xlsx');
    ff_data = xlsread(fn.name{n});
    [W, intmesh, airfoil] = load_intmesh();
    intmesh.thrust_loss = 1;
    intmesh.power_loss = 1;
    intmesh.download = .95;
    
    j = 1;
    for r = ff_data'
        T = r(1)/2;
        V = r(2);
        [~, P, ~] = thrust_bemt(T, 1e-3, intmesh, h(n), V, airfoil, inf);
        P_req(j) = 2*P;
        j=j+1;
    end
    
    V_kts = ff_data(:,2)*0.592484;
    h1(n) = plot(V_kts, P_req)
    hold on
    
    P_avail = 300*density(h(n))/density(0);
    V_max = interp1(P_req, V_kts, P_avail);
    
    yl = ylim;
    y = [P_avail P_avail];
    plot(xlim, y, '--')
    h2(n) = plot(V_max, P_avail, 'o')
    h3(n) = plot([V_max V_max], [yl(1) P_avail], '--')
end


xlabel('V_{\infty} (kts)')
%ylabel('Power Required @ SL (HP)')
ylabel('Power Required (HP)')

legend([h1, h2, h3], 'MCP @ SL', 'MCP @ 3000m', 'V Max @ SL', ...
    'V Max @ 3000m', 'P_{avai} @ SL', 'P_{avai} @ 3000m', 'Location', 'NW')
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
%fprintf('\nMax V_{\infty}: %.2f kts\n', V_max);















