function [lla] = ned2lla(ned, lla0)
    lla=[];
    N = ned(:, 1); 
    E = ned(:, 2);
    D = ned(:, 3);

    p = vpa(pi,20);
    a = 6378137.0; % m
    b = 6356752.3; % m
    ep = sqrt(1 - (b^2 / a^2));

    phi = lla0(2);

    Rm = a * (1 + ep^2 * (2/3 * (sin(phi))^2 - 1));
    Rn = a * (1 + ep^2 / 2 * (sin(phi))^2);

    dmu = N * atan2(1,Rm);
    dl = E * atan2(1,Rn*cos(lla0(1)));
    lla(:, 1) = mod((lla0(1) + dmu)+ p, 2 * p) - p;
    lla(:, 2) = mod((lla0(2) + dl)+ p, 2 * p) - p;
    lla(:, 3) = lla0(3) - D;
end
    