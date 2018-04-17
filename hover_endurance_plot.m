% Run after hover_endurance_sim.m
clc

plot(t_ref/60, P_ref) 
xlabel('Time (min)')
ylabel('Power Required (HP)')
title('Power Required During Hover Endurance Test @ SL')
fprintf('Intermesh Hover Time @ 50%% Fuel:   %d hr %d min\n', h, m);