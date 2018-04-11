function rho = density(h)
    % Calculate Desnity
    % h: feet
    % return: rho in slugs/ft3
    rho = .00238*(1-.00198*h/288.16)^4.2553;
end
