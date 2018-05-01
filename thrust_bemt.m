function [P, P_sl, FM, dT, ru] = thrust_bemt(T_req, epsilon, rotor, h, V_inf, airfoil, max_thr)

    thr = rotor.t;
    T = 0;
    
    while abs(T - T_req) > epsilon && thr < max_thr
        [T, P, P_sl, FM, dT, ru] = bemt(rotor, thr, h, V_inf, airfoil);
        thr = thr + 1/(1+exp((T-T_req)/T_req^(1/1.3)))-.5;
    end
    
    if thr > max_thr
        P = NaN;
        P_sl = NaN;
    end