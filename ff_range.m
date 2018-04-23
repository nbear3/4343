function [R, V_range, P_range] = ff_range(V_kts, P_req)
    y = P_req./V_kts';
    [~, i] = min(y);
    fun = @(x) spline(V_kts, y, x);
    V_range = fminsearch(fun, V_kts(i)); % km/hr

    P_range = interp1(V_kts, P_req, V_range);

    % Parameters
    GW = 600*2.2046;
    phi = 0.55; % includes payload
    sfc = 0.578;
    payload = 100*2.2046;

    % Fuel Fraction Available
    W_fuel = GW-phi*GW-payload; 
    R = W_fuel/(sfc*P_range/V_range)*1.852;    
end