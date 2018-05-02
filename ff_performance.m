clc
clear
close all

MCP = 280;
Dash_THP = MCP*1.1;

h = [0, 3000];

for n = 1:length(h)
    h_ft = h(n)*3.28084; % convert to ft
    [Vkts, P_req] = compound_power(h_ft);
    fprintf('\nPerformance @ %d m\n', h(n));
    
    %% Plot The Max Speed
    figure
    hold on
    P_avail = MCP*density(h_ft)/density(0);
    V_max = interp1(P_req, Vkts, P_avail);
    V_dash = interp1(P_req, Vkts, Dash_THP);
    fprintf('\tMax Speed @ MCP: %.2f kts\n', V_max);
    fprintf('\t@P: %.2f HP\n', P_avail);
    fprintf('\n\tDash Speed: %.2f kts\n', V_dash);
    fprintf('\t@P: %.2f HP\n', Dash_THP);

    yl = ylim;
    plot(Vkts, P_req);
    plot(xlim, [P_avail P_avail], '--')
    plot(V_max, P_avail, 'o');

    
    %% Plot The Max Range
    [R, V_range, P_range] = ff_range(Vkts, P_req);
    fprintf('\n\tTotal Range: %.2f km\n', R);
    fprintf('\t50%% Range: %.2f km\n', R/2);
    fprintf('\t@V: %.2f kts\n', V_range);
    fprintf('\t@P: %.2f HP\n', P_range);
    plot(V_range, P_range, 'o');
    plot([V_max V_max], [yl(1) P_avail], '--');
    plot([0 V_range]*1.2, [0 P_range]*1.2, '--')
    plot([V_range V_range], [0 P_range], '--');
    
    %% Labels
    xlabel('Speed (kts)')
    ylabel('Power Required (HP)')
    legend('P Req', 'MCP @ SL', 'V Max', 'V Best Range', 'Location', 'NW')
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
end

%% Titles

figure(1)
title("Speed and Range @ SL")

figure(2)
title("Speed and Range @ 3000m")













