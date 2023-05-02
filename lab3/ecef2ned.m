function [res] = ecef2ned(ecef, lla)
    % Conversion between the ECEF and NED systems
    % input: ecef (3*N)
    % input: lla (3*N)
    % output: (3*N)

    clat = cos(lla(1));
    slat = sin(lla(1));
    clon = cos(lla(2));
    slon = sin(lla(2));

    R = [-clon*slat, -slon*slat, clat;
         -slon, clon, 0;
         -clon*clat, -slon*clat, -slat];

    res = R * ecef';
end