function [lla] = ecef2lla(inp)
    x = inp(1);
    y = inp(2);
    z = inp(3);

    % WGS84 ellipsoid constants:
    a = 6378137.0;
    e = 8.1819190842622e-2;
    
    % calculations:
    b   = sqrt(a^2*(1-e^2));
    ep  = sqrt((a^2-b^2)/b^2);
    p   = sqrt(x.^2+y.^2);
    th  = atan2(a*z,b*p);
    lon = atan2(y,x);
    lat = atan2((z+ep^2.*b.*sin(th).^3),(p-e^2.*a.*cos(th).^3));
    N   = a./sqrt(1-e^2.*sin(lat).^2);
    alt = p./cos(lat)-N;
    
    lon = mod(lon,2*pi);
    
    k=abs(x)<1 & abs(y)<1;
    alt(k) = abs(z(k))-b;
    lla=[rad2deg(lat) rad2deg(lon) alt];
end