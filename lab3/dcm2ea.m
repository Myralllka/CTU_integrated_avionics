function [res] = dcm2ea(dcm)
    % Conversion of DCM to  Euler angles
    phi = atan2(dcm(3, 2), dcm(3, 3));
    the = -asin(dcm(3, 1));
    psi = atan2(dcm(2, 1), dcm(1, 1));

    res = [phi the psi];

end