clear all;

% Data fusion using Kalman filtering
% 1) Implement Kalman filter for position estimation, velocity
% and orientation (euler angles) in Matlab. The model for 
% implementation is described in the file "Nav-assignment-v1-2023.pdf" on Moodle.
% 
% 2) Kalman filtering results will be represented graphically as:
%   a) Position in geographic system LLA from KF and GNSS;
%   b) Velocity in NED and corresponding velocity from GNSS;
%   c) Euler position angles (pitch roll yaw) obtained from:
%      i)   Three-axis accelerometer;
%      ii)  Gyroscopes using dara integration corrected for estimated bias;
%      iii) KF-state vector;
%   d) Kalman gain (of all parts of state vector);
%   e) Innovation and their histogram
%   f) P Matrices (diagonal elements)

% Data description:
% Dataset: 2022-05-13_CarData.mat
% Sampling frequency - IMU - 100 Hz
% Sampling frequency - GNSS - 10 Hz
% Initial heading: 177deg.
% Data structure:
% Inertial data:
% RawIMU.SF:
% 01 - accelerometer data X – (g) – longitudinal
% 02 - accelerometer data Y – (g) – lateral
% 03 - accelerometer data Z – (g) – vertical
% RawIMU.W:
% 01 - gyroscope data X – (rad) - longitudinal
% 02 - gyroscope data Y – (rad) – lateral
% 03 - gyroscope data Z – (rad) – vertical
% RawGNSS.LLA:
% 01 - Latitude – (°)
% 02 - Longitude – (°)
% 03 - Altitude – (m)
% RawGNSS.VEL:
% 01 – GNSS Velocity X (ECEF) (m/s)
% 02 – GNSS Velocity Y (ECEF) (m/s)
% 03 – GNSS Velocity Z (ECEF) (m/s)

data = load('2022-05-13_CarData.mat');
dt = 1/100;          % IMU freq 100Hz
dt_gnss = 1/10;      % GNSS freq 10Hz
%%

% Gyro data normalization
mean_gx = mean(data.RawIMU.W(1:1000, 1));
mean_gy = mean(data.RawIMU.W(1:1000, 2));
mean_gz = mean(data.RawIMU.W(1:1000, 3));
data.RawIMU.W(:, 1) = data.RawIMU.W(:, 1) - mean_gx;
data.RawIMU.W(:, 2) = data.RawIMU.W(:, 2) - mean_gy;
data.RawIMU.W(:, 3) = data.RawIMU.W(:, 3) - mean_gz;

% IMU acceleration pre-processing
% g_multiplier = mean(sqrt(sum(data.RawIMU.SF(1:1000, :).^2, 2)));
g = comp_gravity(data.RawGNSS.LLA(2, :), "deg");
% g_multiplier = g_multiplier / g;
mean_ax = mean(data.RawIMU.SF(1:1000, 1) * g);
mean_ay = mean(data.RawIMU.SF(1:1000, 2) * g);
mean_az = mean(data.RawIMU.SF(1:1000, 3) * g);

roll_init_acc = mean(atan(data.RawIMU.SF(1:1000, 2) ./ data.RawIMU.SF(1:1000, 3))); % acc angles
pitch_init_acc = mean(-atan(data.RawIMU.SF(1:1000, 1) ./ data.RawIMU.SF(1:1000, 3)));

init_angle = [roll_init_acc, pitch_init_acc, deg2rad(177)];
Rbn = ea2dcm(init_angle);

mech_position = [[0, 0, 0]];
mech_velocity = [[0, 0, 0]];
a_centripetal = [[0, 0, 0]];

angles_acc = [[roll_init_acc, pitch_init_acc]];
angles_gyro = [[roll_init_acc, pitch_init_acc, init_angle(3)]];
angles_gnss = [[init_angle(3)]];

mech_angles = [[0, 0, 0]];

lla0 = deg2rad(data.RawGNSS.LLA(2, :));
gps_track = [];
N = size(data.RawGNSS.LLA, 1);

%% Kalman filter: data fusion
kalman_position = [[0, 0, 0]];

O = [0 0 0; 0 0 0; 0 0 0];
I = eye(3);

% Values initialisation
Fpv = eye(3);
Fpv(3, 3) = -1;
Fpp = [0 0 0; 0 0 0; 0 0 0];
Fpr = [0 0 0; 0 0 0; 0 0 0];
Fvp = [0 0 0; 0 0 0; 0 0 0];
Fvv = [0 0 0; 0 0 0; 0 0 0];
Fvr = Rbn * sq(data.RawIMU.SF(1, :).' * g);
Frp = [0 0 0; 0 0 0; 0 0 0];
Frv = [0 0 0; 0 0 0; 0 0 0];
Frr = [0 0 0; 0 0 0; 0 0 0];

A = [Fpp, Fpv, Fpr, O, O;
     Fvp, Fvv, Fvr, -Rbn, O;
     Frp, Frv, Frr, O, Rbn;
     O, O, O, O, O;
     O, O, O, O, O];
G = [I O O O O; 
     O -Rbn O O O; 
     O O Rbn O O;
     O O O I O;
     O O O O I];
I15 = eye(15);
xk = zeros(1, 15);
Pk = eye(15);

Qc = eye(15);

Qk = G * Qc * G.' * dt;
Rk = diag([1.5, 1.5, 1.5, 0.1, 0.1, 0.1])^2;
% Rk = diag([0.1 0.1 0.1 0.1 0.1 0.1])^2;
% Rk = diag([0 0 0 0 0 0])^2;
H = [I O O O O; 
     O I O O O];
Qc = eye(15);
Fi = I15 + A*dt;

inn_array = zeros(1, 6);
pss_array = [Pk];
kss_array = [];

% Kalman filter:
for i = 1001:N;
    Fvr = Rbn * sq(data.RawIMU.SF(1, :).' * g);
    A = [Fpp, Fpv, Fpr, O, O;
         Fvp, Fvv, Fvr, -Rbn, O;
         Frp, Frv, Frr, O, Rbn;
         O, O, O, O, O;
         O, O, O, O, O];
    Fi = I15 + A*dt;

    G = [I O O O O; 
         O -Rbn O O O; 
         O O Rbn O O;
         O O O I O;
         O O O O I];
    Qk = G * Qc * G.' * dt;

    % => Time update, predict
    % -> Extrapolate the state
    x_next = Fi * xk(end, :).';
    % -> Extrapolate the uncertainty
    P_next = Fi * Pk * Fi.' + Qk;

    % => Mechanisation
    g_tmp = comp_gravity(data.RawGNSS.LLA(i, :), "deg");
    if ~isnan(g_tmp)
        g = g_tmp;
    end
    omega = data.RawIMU.W(i, :);
    ac = a_centripetal(end, :).';
    vp = (Rbn * (data.RawIMU.SF(i, :).' * g - ac) + [0; 0; g]).';
    pp = mech_velocity(end, :);
    pp(3) = -pp(3);
    Rbnp = Rbn * sq(omega);
    vb = Rbn.' * mech_velocity(end, :).';
    a_centripetal(end+1, :) = (sq(omega) * vb).';
    mech_position(end+1, :) = mech_position(end, :) + pp*dt;
    mech_velocity(end+1, :) = mech_velocity(end, :) + vp*dt;
    Rbn = norm_DCM(Rbn + Rbnp*dt);
    mech_angles(end+1, :) = dcm2ea(Rbnp);
    phi = angles_gyro(end, 1);
    theta = angles_gyro(end, 2);
    R = [1,  sin(phi)*tan(theta),    cos(phi) * tan(theta);
         0,  cos(phi),               -sin(phi);
         0,  sin(phi)/cos(theta),    cos(phi)/cos(theta)];
    om = data.RawIMU.W(i, :);
    rpy_p = (R * om.').';
    angles_gyro(end+1, :) = angles_gyro(end, :) + rpy_p*dt;
    
    % -> zk compute
    angls = [deg2rad(data.RawGNSS.LLA(i, 1)), deg2rad(data.RawGNSS.LLA(i, 2)), data.RawGNSS.LLA(i, 3)];
    ned_gnss_vel = ecef2ned(data.RawGNSS.VEL(i, :), angls);
    if ~isnan(atan2(ned_gnss_vel(2), ned_gnss_vel(1)))
        gps_track(end+1, :) = lla2ned(deg2rad(data.RawGNSS.LLA(i, :)).', lla0.');
        % => Kompute the Kalman gain
        K = P_next * H.' * inv(H * P_next * H.'  + Rk);
        % 
        p_gnss = lla2ned(deg2rad(data.RawGNSS.LLA(i, :)).',lla0.');
        v_gnss = ned_gnss_vel;
        p_mech = mech_position(end, :).';
        v_mech = mech_velocity(end, :).';
        zk = [p_gnss - p_mech; v_gnss - v_mech];
        expected = H * x_next;
        inn_array(end + 1, :) = (zk - expected).';
        xk(end + 1, :) = (x_next + K * inn_array(end, :).').';

        % Remark from the lecture
        % K(K>1) = 1;
        % K(K<-1) = -1;

        % Update position:
        kalman_position(end+1, :) = mech_position(end, :).' + xk(end, 1:3).';

        % P - symetric positive definite!
        Pk = (I15 - K*H) * P_next;
        Pk = (Pk + Pk.') / 2;
        pss_array(i, :) = diag(Pk);
        kss_array(:, :, i) = K;

        % 
        mech_position(end, :) = mech_position(end, :) + xk(end, 1:3);
        mech_velocity(end, :) = mech_velocity(end, :) + xk(end, 4:6);
        Rbn = (I - ea2dcm(xk(end, 7:9)));
        xk(end, :) = 0;

    end

end
% x-, Pk- - extrapolation
% x+, Pk+ - observation
p_gnss = [0 0 0];
v_gnss = [0 0 0];
for i = 1:N
end

figure(1);
plot(gps_track(:, 2), gps_track(:, 1));
hold on
% plot(mech_position(:, 2), mech_position(:, 1));
% grid on
plot(kalman_position(:, 2), kalman_position(:, 1));
grid on
legend("GPS", "Kalman filter")
