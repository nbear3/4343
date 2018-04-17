clc
clear 

% rotor parameter
rotor.OR = 650;
rotor.t = 5;
rotor.R = 1.3*3.28084;
rotor.c = 0.8;
rotor.tr = 0.47;
rotor.n = 2;
rotor.thrust_loss = cosd(13);
rotor.power_loss = 1.1;
rotor.download = 0.9;

% naca23012
airfoil.root = [-10.7500000000000,-0.777600000000000,0.0553800000000000;-10.5000000000000,-0.789400000000000,0.0523100000000000;-10.2500000000000,-0.800500000000000,0.0483300000000000;-10,-0.813700000000000,0.0432200000000000;-9.75000000000000,-0.831900000000000,0.0364400000000000;-9.50000000000000,-0.838900000000000,0.0313800000000000;-9.25000000000000,-0.833300000000000,0.0284900000000000;-9,-0.821300000000000,0.0260700000000000;-8.75000000000000,-0.804100000000000,0.0250300000000000;-8.50000000000000,-0.786000000000000,0.0242400000000000;-8.25000000000000,-0.768500000000000,0.0233000000000000;-8,-0.751300000000000,0.0223200000000000;-7.75000000000000,-0.734000000000000,0.0212800000000000;-7.50000000000000,-0.716400000000000,0.0203000000000000;-7.25000000000000,-0.698200000000000,0.0195700000000000;-7,-0.680200000000000,0.0180500000000000;-6.75000000000000,-0.661400000000000,0.0174100000000000;-6.50000000000000,-0.642000000000000,0.0168900000000000;-6.25000000000000,-0.605500000000000,0.0163000000000000;-6,-0.568300000000000,0.0157000000000000;-5.75000000000000,-0.532700000000000,0.0148000000000000;-5.50000000000000,-0.497400000000000,0.0140800000000000;-5.25000000000000,-0.459700000000000,0.0135900000000000;-5,-0.420400000000000,0.0131800000000000;-4.75000000000000,-0.387300000000000,0.0124700000000000;-4.50000000000000,-0.349100000000000,0.0120300000000000;-4.25000000000000,-0.311200000000000,0.0116800000000000;-4,-0.281000000000000,0.0111200000000000;-3.75000000000000,-0.249500000000000,0.0108500000000000;-3.50000000000000,-0.221800000000000,0.0103800000000000;-3.25000000000000,-0.194600000000000,0.0101400000000000;-3,-0.168800000000000,0.00982000000000000;-2.75000000000000,-0.143100000000000,0.00956000000000000;-2.50000000000000,-0.117500000000000,0.00930000000000000;-2.25000000000000,-0.0970000000000000,0.00829000000000000;-2,-0.0785000000000000,0.00727000000000000;-1.75000000000000,-0.0559000000000000,0.00688000000000000;-1.50000000000000,-0.0320000000000000,0.00667000000000000;-1.25000000000000,-0.00770000000000000,0.00654000000000000;-1,0.0170000000000000,0.00648000000000000;-0.750000000000000,0.0419000000000000,0.00647000000000000;-0.500000000000000,0.0669000000000000,0.00650000000000000;-0.250000000000000,0.0923000000000000,0.00657000000000000;0,0.117500000000000,0.00668000000000000;0.250000000000000,0.142500000000000,0.00683000000000000;0.500000000000000,0.168500000000000,0.00703000000000000;0.750000000000000,0.196000000000000,0.00723000000000000;1,0.225000000000000,0.00743000000000000;1.25000000000000,0.257600000000000,0.00765000000000000;1.50000000000000,0.290100000000000,0.00789000000000000;1.75000000000000,0.325600000000000,0.00814000000000000;2,0.364400000000000,0.00840000000000000;2.25000000000000,0.400700000000000,0.00864000000000000;2.50000000000000,0.437200000000000,0.00889000000000000;2.75000000000000,0.475800000000000,0.00912000000000000;3,0.513000000000000,0.00934000000000000;3.25000000000000,0.549600000000000,0.00956000000000000;3.50000000000000,0.588000000000000,0.00974000000000000;3.75000000000000,0.626000000000000,0.00991000000000000;4.25000000000000,0.674100000000000,0.0102100000000000;4.50000000000000,0.696700000000000,0.0103900000000000;4.75000000000000,0.719400000000000,0.0105700000000000;5,0.742300000000000,0.0107100000000000;5.25000000000000,0.764400000000000,0.0109500000000000;5.50000000000000,0.787200000000000,0.0111200000000000;5.75000000000000,0.810100000000000,0.0112800000000000;6,0.832000000000000,0.0115500000000000;6.25000000000000,0.854900000000000,0.0117400000000000;6.50000000000000,0.878100000000000,0.0119100000000000;6.75000000000000,0.900400000000000,0.0121800000000000;7,0.923400000000000,0.0124000000000000;7.25000000000000,0.947000000000000,0.0125700000000000;7.50000000000000,0.969200000000000,0.0128700000000000;7.75000000000000,0.992700000000000,0.0130800000000000;8,1.01630000000000,0.0132900000000000;8.25000000000000,1.03860000000000,0.0136400000000000;8.50000000000000,1.06320000000000,0.0138000000000000;8.75000000000000,1.08610000000000,0.0141300000000000;9,1.10990000000000,0.0143900000000000;9.25000000000000,1.13320000000000,0.0147000000000000;9.50000000000000,1.15630000000000,0.0150400000000000;9.75000000000000,1.17930000000000,0.0153900000000000;10,1.20210000000000,0.0157600000000000;10.2500000000000,1.22350000000000,0.0162400000000000;10.5000000000000,1.24610000000000,0.0166200000000000;10.7500000000000,1.26680000000000,0.0171500000000000;11,1.28720000000000,0.0176900000000000;11.2500000000000,1.30720000000000,0.0182500000000000;11.5000000000000,1.32550000000000,0.0189200000000000;11.7500000000000,1.34280000000000,0.0196400000000000;12,1.36000000000000,0.0203400000000000;12.2500000000000,1.37420000000000,0.0212100000000000;12.5000000000000,1.38700000000000,0.0220300000000000;12.7500000000000,1.39670000000000,0.0229600000000000;13,1.40540000000000,0.0239900000000000;13.2500000000000,1.41240000000000,0.0251800000000000;13.5000000000000,1.41980000000000,0.0264100000000000;13.7500000000000,1.42060000000000,0.0282000000000000;14,1.42790000000000,0.0296400000000000;14.2500000000000,1.42990000000000,0.0316200000000000;14.5000000000000,1.42410000000000,0.0345500000000000;14.7500000000000,1.42570000000000,0.0370100000000000;15,1.42420000000000,0.0400100000000000;15.2500000000000,1.41900000000000,0.0436900000000000;15.5000000000000,1.40870000000000,0.0482900000000000;15.7500000000000,1.39320000000000,0.0538300000000000;16,1.37300000000000,0.0602200000000000;16.2500000000000,1.35610000000000,0.0662500000000000;16.5000000000000,1.33940000000000,0.0723100000000000;16.7500000000000,1.32110000000000,0.0786600000000000;17,1.30280000000000,0.0851200000000000];
airfoil.tip = [-15.2500000000000,-1.06270000000000,0.0622000000000000;-15,-1.09520000000000,0.0531300000000000;-14.7500000000000,-1.11590000000000,0.0470600000000000;-14.5000000000000,-1.12910000000000,0.0427000000000000;-14.2500000000000,-1.13810000000000,0.0393200000000000;-14,-1.14260000000000,0.0367200000000000;-13.7500000000000,-1.14280000000000,0.0347400000000000;-13.5000000000000,-1.13960000000000,0.0332000000000000;-13.2500000000000,-1.13290000000000,0.0321100000000000;-13,-1.14690000000000,0.0289000000000000;-12.7500000000000,-1.13580000000000,0.0277000000000000;-12.5000000000000,-1.12010000000000,0.0269100000000000;-12.2500000000000,-1.10350000000000,0.0261400000000000;-12,-1.08620000000000,0.0253400000000000;-11.7500000000000,-1.06870000000000,0.0245000000000000;-11.5000000000000,-1.05070000000000,0.0236400000000000;-11.2500000000000,-1.03230000000000,0.0227800000000000;-11,-1.01330000000000,0.0219500000000000;-10.7500000000000,-0.993800000000000,0.0211900000000000;-10.5000000000000,-0.973700000000000,0.0205600000000000;-10.2500000000000,-0.956300000000000,0.0194800000000000;-10,-0.939400000000000,0.0184500000000000;-9.75000000000000,-0.919600000000000,0.0179700000000000;-9.50000000000000,-0.899600000000000,0.0175400000000000;-9.25000000000000,-0.880100000000000,0.0170600000000000;-9,-0.861000000000000,0.0165700000000000;-8.75000000000000,-0.842200000000000,0.0160800000000000;-8.50000000000000,-0.823500000000000,0.0156300000000000;-8.25000000000000,-0.804600000000000,0.0152300000000000;-8,-0.783800000000000,0.0149200000000000;-7.75000000000000,-0.753200000000000,0.0138500000000000;-7.50000000000000,-0.719100000000000,0.0134400000000000;-7.25000000000000,-0.684700000000000,0.0130300000000000;-7,-0.650200000000000,0.0126100000000000;-6.75000000000000,-0.617300000000000,0.0122300000000000;-6.50000000000000,-0.583300000000000,0.0117800000000000;-6.25000000000000,-0.549100000000000,0.0112000000000000;-6,-0.513000000000000,0.0108800000000000;-5.75000000000000,-0.481800000000000,0.0105800000000000;-5.50000000000000,-0.449800000000000,0.0103100000000000;-5.25000000000000,-0.421200000000000,0.00985000000000000;-5,-0.392500000000000,0.00960000000000000;-4.75000000000000,-0.366000000000000,0.00940000000000000;-4.50000000000000,-0.340400000000000,0.00914000000000000;-4.25000000000000,-0.315100000000000,0.00894000000000000;-4,-0.289200000000000,0.00883000000000000;-3.75000000000000,-0.263200000000000,0.00862000000000000;-3.50000000000000,-0.236400000000000,0.00853000000000000;-3.25000000000000,-0.210000000000000,0.00838000000000000;-3,-0.183100000000000,0.00832000000000000;-2.75000000000000,-0.156400000000000,0.00821000000000000;-2.50000000000000,-0.129600000000000,0.00814000000000000;-2.25000000000000,-0.102600000000000,0.00805000000000000;-2,-0.0761000000000000,0.00792000000000000;-1.75000000000000,-0.0554000000000000,0.00692000000000000;-1.50000000000000,-0.0316000000000000,0.00644000000000000;-1.25000000000000,-0.00630000000000000,0.00619000000000000;-1,0.0194000000000000,0.00603000000000000;-0.750000000000000,0.0454000000000000,0.00595000000000000;-0.500000000000000,0.0716000000000000,0.00591000000000000;-0.250000000000000,0.0977000000000000,0.00587000000000000;0,0.124100000000000,0.00587000000000000;0.250000000000000,0.150300000000000,0.00591000000000000;0.500000000000000,0.176700000000000,0.00597000000000000;0.750000000000000,0.202500000000000,0.00603000000000000;1,0.228000000000000,0.00612000000000000;1.25000000000000,0.254400000000000,0.00627000000000000;1.50000000000000,0.281000000000000,0.00642000000000000;1.75000000000000,0.307600000000000,0.00657000000000000;2,0.334800000000000,0.00674000000000000;2.25000000000000,0.361700000000000,0.00691000000000000;2.50000000000000,0.390900000000000,0.00709000000000000;2.75000000000000,0.420200000000000,0.00724000000000000;3,0.453100000000000,0.00742000000000000;3.25000000000000,0.484700000000000,0.00761000000000000;3.50000000000000,0.518800000000000,0.00775000000000000;3.75000000000000,0.555300000000000,0.00794000000000000;4,0.591300000000000,0.00808000000000000;4.25000000000000,0.626100000000000,0.00828000000000000;4.50000000000000,0.660300000000000,0.00840000000000000;4.75000000000000,0.694900000000000,0.00857000000000000;5,0.730600000000000,0.00873000000000000;5.25000000000000,0.766000000000000,0.00886000000000000;5.50000000000000,0.800800000000000,0.00905000000000000;5.75000000000000,0.831300000000000,0.00918000000000000;6,0.854500000000000,0.00930000000000000;6.25000000000000,0.877100000000000,0.00948000000000000;6.50000000000000,0.900000000000000,0.00963000000000000;6.75000000000000,0.923200000000000,0.00977000000000000;7,0.945800000000000,0.00997000000000000;7.25000000000000,0.968900000000000,0.0101500000000000;7.50000000000000,0.992200000000000,0.0103100000000000;7.75000000000000,1.01460000000000,0.0105800000000000;8,1.03850000000000,0.0107300000000000;8.25000000000000,1.06160000000000,0.0109700000000000;8.50000000000000,1.08480000000000,0.0112100000000000;8.75000000000000,1.10820000000000,0.0114500000000000;9,1.13130000000000,0.0117400000000000;9.25000000000000,1.15450000000000,0.0120600000000000;9.50000000000000,1.17800000000000,0.0123600000000000;9.75000000000000,1.20090000000000,0.0127400000000000;10,1.22430000000000,0.0130700000000000;10.2500000000000,1.24680000000000,0.0135000000000000;10.5000000000000,1.26890000000000,0.0139600000000000;10.7500000000000,1.29130000000000,0.0143800000000000;11,1.31280000000000,0.0148900000000000;11.2500000000000,1.33360000000000,0.0154400000000000;11.5000000000000,1.35470000000000,0.0159400000000000;11.7500000000000,1.37460000000000,0.0165300000000000;12.2500000000000,1.41280000000000,0.0177700000000000;12.5000000000000,1.43050000000000,0.0184400000000000;12.7500000000000,1.44750000000000,0.0191300000000000;13,1.46360000000000,0.0198300000000000;13.2500000000000,1.47660000000000,0.0205800000000000;13.5000000000000,1.48750000000000,0.0214000000000000;13.7500000000000,1.49870000000000,0.0222200000000000;14,1.50580000000000,0.0233400000000000;14.2500000000000,1.51570000000000,0.0243400000000000;14.5000000000000,1.52390000000000,0.0255000000000000;14.7500000000000,1.52870000000000,0.0269900000000000;15,1.53040000000000,0.0288600000000000;15.2500000000000,1.53680000000000,0.0304800000000000;15.5000000000000,1.54040000000000,0.0324800000000000;15.7500000000000,1.54130000000000,0.0349000000000000;16,1.53900000000000,0.0378900000000000;16.2500000000000,1.53250000000000,0.0416000000000000;16.5000000000000,1.52020000000000,0.0463600000000000;16.7500000000000,1.50260000000000,0.0521200000000000;17,1.48570000000000,0.0580300000000000;17.2500000000000,1.46490000000000,0.0646000000000000;17.5000000000000,1.43780000000000,0.0722700000000000;17.7500000000000,1.40510000000000,0.0809100000000000;18,1.37020000000000,0.0900200000000000];
airfoil.ratio = .5;

% Parameters
GW = 600*2.2046;
phi = 0.55; % includes payload
sfc = 0.578;

% Mission Parameters
payload = 100*2.2046; % kg

% Performance Requirements

h = 0*3000*3.28084; % height (ft)
V = 0; % forward speed (ft/s)

% Fuel Fraction Available
Rf_avail = 1-phi-payload./GW; 

W = GW;    
Rf_50 = Rf_avail/2;
Rf_next = 0;
epsilon = 1e-3;
max_thr = 30; % maximum collective pitch

t = 0;
time_step = 100; % in seconds

i = 1;
while (Rf_avail - Rf_next) > Rf_50
    [~, P_req, ~] = thrust_bemt(W/2, epsilon, rotor, h, V, airfoil, max_thr);
    F_req = 2*P_req*time_step*sfc/3600; % requires 2*P_req for both rotor
    W = W - F_req;

    Rf_next = F_req/GW;
    Rf_avail = Rf_avail - Rf_next;
    t = t + time_step;  
    
    F_ref(i) = F_req;
    P_ref(i) = 2*P_req;
    t_ref(i) = t;
    i = i + 1;
end

F_left = (Rf_avail - Rf_50)*GW;
[~, P_req, ~] = thrust_bemt(W/2, epsilon, rotor, h, V, airfoil, max_thr);
t = t + F_left/(2*P_req)/sfc*3600;
h = floor(t/3600);
m = floor(mod(t, 3600)/60);
fprintf('Intermesh Hover Time @ 50%% Fuel:   %d hr %d min\n', h, m);