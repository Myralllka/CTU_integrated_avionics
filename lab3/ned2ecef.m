function ecef = ned2ecef(ned, lla0)
    % Conversion between the NED and ECEF systems
    slat = sin(lla0(1));
    clat = cos(lla0(1));
    slon = sin(lla0(2));
    clon = cos(lla0(2));
    
    R = [  -slat*clon  -slat*slon   clat;
           -slon          clon         0;
           -clat*clon  -clat*slon  -slat];
    difned = R' * ned';

    ecef = difned';
    % ecef = (R' * ned')';
end
