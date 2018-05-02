clc
clear
close all

h = [0, 3000];

for n = 1:length(h)
    h_ft = h(n)*3.28084; % convert to ft
    [Vkts, P_req] = compound_power(h_ft);
    
%     figure(n)
    h1(n) = plot(Vkts, P_req);
    hold on

    P_avail = 280*density(h_ft)/density(0);
    V_max = interp1(P_req, Vkts, P_avail);
    fprintf('\nMaximum MCP Speed @ %d m\n', h(n));
    fprintf('\tMax Speed: %.2f kts\n', V_max);
    fprintf('\t@P: %.2f HP\n', P_avail);
    
    [R, V_range, P_range] = ff_range(Vkts, P_req);
    fprintf('\nBest Range @ %d m\n', h(n));
    fprintf('\tTotal Range: %.2f km\n', R);
    fprintf('\t50%% Range: %.2f km\n', R/2);
    fprintf('\t@V: %.2f kts\n', V_range);
    fprintf('\t@P: %.2f HP\n', P_range);

    yl = ylim;
    y = [P_avail P_avail];
    plot(xlim, y, '--')
    h2(n) = plot(V_max, P_avail, 'o');
    h3(n) = plot([V_max V_max], [yl(1) P_avail], '--');
    h4(n) = plot(V_range, P_range, 'x');
end


xlabel('Speed (kts)')
ylabel('Power Required (HP)')
legend([h1, h2, h3], 'P Req @ SL', 'P Req @ 3000m', 'V Max @ SL', ...
    'V Max @ 3000m', 'MCP @ SL', 'MCP @ 3000m', 'Location', 'NW')
set(findall(gca, 'Type', 'Line'),'LineWidth',1);















