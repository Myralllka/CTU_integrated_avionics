function [res] = ea2dcm(ea)
    % Conversion of Euler angles to DCM
    % Roll Pitch Yaw to DCM
    
    cphi = cos(ea(1));
    sphi = sin(ea(1));
    cthe = cos(ea(2));
    sthe = sin(ea(2));
    cpsi = cos(ea(3));
    spsi = sin(ea(3));
    res = [ cthe*cpsi, -cphi*spsi + sphi*sthe*cpsi, sphi*spsi + cphi*sthe*cpsi;
            cthe*spsi, cphi*cpsi + sphi*sthe*spsi, -sphi*cpsi + cphi*sthe*spsi;
            -sthe, sphi*cthe, cphi*cthe];
    % a = ea(1); % roll
    % b = ea(2); % pitch
    % c = ea(3); % yaw

    % Ax = [1 0 0;
    %     0 cos(a) sin(a);
    %     0 -sin(a) cos(a)];
    % Ay = [cos(b) 0 -sin(b);
    %     0 1 0;
    %     sin(b) 0 cos(b)];
    % Az = [cos(c) sin(c) 0;
    %     -sin(c) cos(c) 0;
    %     0 0 1];
    % res = (Az * Ay * Ax)';

end
