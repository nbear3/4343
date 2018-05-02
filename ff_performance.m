clc
clear
close all

MCP = 280;
h = [0, 3000];

for n = 1:length(h)
    h_ft = h(n)*3.28084; % convert to ft
    [Vkts, P_req] = compound_power(h_ft);
    
    %% Plot The Max Speed
    figure
    hold on
    P_avail = MCP*density(h_ft)/density(0);
    V_max = interp1(P_req, Vkts, P_avail);
    fprintf('\nMaximum MCP Speed @ %d m\n', h(n));
    fprintf('\tMax Speed: %.2f kts\n', V_max);
    fprintf('\t@P: %.2f HP\n', P_avail);

    yl = ylim;
    plot(Vkts, P_req);
    plot(xlim, [P_avail P_avail], '--')
    plot(V_max, P_avail, 'o');
    plot([V_max V_max], [yl(1) P_avail], '--');
    xlabel('Speed (kts)')
    ylabel('Power Required (HP)')
    legend('P Req', 'MCP @ SL', 'V Max', 'Location', 'NW')
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
    
    %% Plot The Max Range
    figure
    hold on
    [R, V_range, P_range] = ff_range(Vkts, P_req);
    fprintf('\nBest Range @ %d m\n', h(n));
    fprintf('\tTotal Range: %.2f km\n', R);
    fprintf('\t50%% Range: %.2f km\n', R/2);
    fprintf('\t@V: %.2f kts\n', V_range);
    fprintf('\t@P: %.2f HP\n', P_range);
    plot(Vkts, P_req);
    plot(V_range, P_range, 'o');
    plot([0 V_range]*1.2, [0 P_range]*1.2, '--')
    plot([V_range V_range], [0 P_range], '--');
    legend('P Req', 'V Best Range', 'Location', 'NW')
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
end

%% Titles

figure(1)
title("Max Speed @ SL")

figure(2)
title("Max Range @ SL")

figure(3)
title("Max Speed @ 3000m")

figure(4)
title("Max Range @ 3000m")












