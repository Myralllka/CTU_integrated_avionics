function [res] = ea2dcm(ea)
    % Conversion of Euler angles to DCM
    % Roll Pitch Yaw to DCM
    
    cphi = cos(ea(1));
    sphi = sin(ea(1));
    cthe = cos(ea(2));
    sthe = sin(ea(2));
    cpsi = cos(ea(3));
    spsi = sin(ea(3));
    
    res = [ cthe*cpsi, -cphi*spsi + sphi*sthe*spsi, sphi*spsi + cphi*sthe*cpsi;
            cthe*spsi, cphi*cpsi + sphi*sthe*spsi, -sphi * cpsi + cphi*sthe*spsi;
            -sthe, sphi * cthe, cphi*cthe];

end
