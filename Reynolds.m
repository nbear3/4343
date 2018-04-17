%p = [0.00238, 0.00176420];
p = 0.00238; % density (slugs/ft^3)
c = 0.79; % chord in ft
R = 4.26; % R in ft
tr = 0.47; % taper ratio

sn = 100;
rp = linspace(1/sn, 1, sn); % Radius Proportion
ru = R * rp; % Radius in Desired Units

cr = 2 * c / (tr + 1); % Chord at Root
ct = cr * tr; % Chord at Tip
c = cr - (cr - ct) * rp; % Chord along Span
u = p*1.5111E-5*3.28084^2;
V = linspace(0,650,sn);
Re = p*V.*c/u;

ru = ru/3.28084; % radius in m
plot(ru, Re)
xlabel('Radius Span (m)')
ylabel('Reynolds Number')