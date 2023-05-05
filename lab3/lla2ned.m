function [res] = lla2ned(lla, lla0)
    % Conversion of LLA to NED
    res = [];
    a = 6378137.0; % m
    b = 6356752.3; % m
    ep = sqrt(1 - (b^2 / a^2));

    phi = lla(1);

    Rm = a * (1 + ep^2 * (2/3 * (sin(phi))^2 - 1));
    Rn = a * (1 + ep^2 / 2 * (sin(phi))^2);

    res(1, 1) = (lla(1, 1) - lla0(1)) / atan(1 / Rm);
    res(2, 1) = (lla(2, 1) - lla0(2)) / atan(1 / (Rn*cos(lla0(1))));
    res(3, 1) = -lla(3, 1) + lla0(3);
end
