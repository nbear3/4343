function rho = density(h)
    % Returns density (rho) at a given altitude (h)
    % rho: desnity in kg/m3
    % h: altitude in m
    % From https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
    
%     rho = 101.29/.2869*(288.14-.00649*h)^4.256/288.08^5.256;
    rho = 2116/1718*(518.7-.00356*h)^4.256/518.6^5.256;