clear all;
close all;

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
gnss_velocity = [[0, 0, 0]];
a_centripetal = [[0, 0, 0]];

angles_acc = [[roll_init_acc, pitch_init_acc]];
angles_gyro = [[roll_init_acc, pitch_init_acc, init_angle(3)]];
angles_gnss = [[init_angle(3)]];

mech_angles = [[0, 0, 0]];

lla0 = deg2rad(data.RawGNSS.LLA(2, :));
gps_track = [];
gps_track_lla = [];
N = size(data.RawGNSS.LLA, 1);

%% Kalman filter: data fusion
kalman_position = [[0, 0, 0]];
kalman_angles = [[0, 0, 0]];
kalman_velocity = [[0, 0, 0]];
kalman_position_lla = [];

O = [0 0 0; 0 0 0; 0 0 0];
I = eye(3);

% Values initialisation
Fpv = eye(3);
% Fpv(3, 3) = -1;
Fpp = [0 0 0; 0 0 0; 0 0 0];
Fpr = [0 0 0; 0 0 0; 0 0 0];
Fvp = [0 0 0; 0 0 0; 0 0 0];
Fvv = [0 0 0; 0 0 0; 0 0 0];
Fvr = -sq(Rbn * (data.RawIMU.SF(1, :).' * g));
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
xk_history = zeros(1, 15);
xk(1, 7:9) = init_angle;
xk(1, 13:15) = mean(data.RawIMU.W(1:1000, :));

Pk = eye(15);
P_next = Pk;
x_next = xk(end, :).';

Qc = ones(15, 1);
Qc(1:3) = 0.5;
Qc(4:6) = 0.05;
Qc(7:9) = 0.05;
Qc(10:12) = 0.0001;
Qc(13:15) = 0.0001;
Qc = diag(Qc);

Rk = diag([1.5, 1.5, 1.5, 0.005, 0.005, 0.005])^2;
H = [I O O O O; 
     O I O O O];
Fi = I15 + A*dt;
Qk = 0.5 * dt * (Fi * G * Qc * G' + G * Qc * G' * Fi');

inn_array = zeros(1, 6);
pss_array = [diag(Pk).'];
kss_array = [];
biases_gyro = [0 0 0];
biases_acc = [0 0 0];
bias_acc = [0 0 0];
bias_gyro = [0 0 0];
K_diag_hist = zeros([6, N]);

% Kalman filter:
for i = 1001:N;
	if mod(i, 10000) == 0;
		i
	end

    % Save values to the history:
    kalman_position(end+1, :) = mech_position(end, :) + xk(end, 1:3);
    kalman_velocity(end+1, :) = mech_velocity(end, :) + xk(end, 4:6);
    kalman_angles(end+1, :) = dcm2ea(Rbn);
    biases_acc(end, :) = bias_acc;
    biases_gyro(end, :) = bias_gyro;
    pss_array(i, :) = diag(Pk);

    % Calculation of pitch and roll from triaxial accelerometer data
    roll = atan(data.RawIMU.SF(i, 2) / data.RawIMU.SF(i, 3)); % acc angles
    pitch = -atan(data.RawIMU.SF(i, 1) / data.RawIMU.SF(i, 3));
    angles_acc(end+1, :) = [roll, pitch];
	
	% IMU data processing
    angls = [deg2rad(data.RawGNSS.LLA(i, 1)), deg2rad(data.RawGNSS.LLA(i, 2)), data.RawGNSS.LLA(i, 3)];
    ned_gnss_vel = ecef2ned(data.RawGNSS.VEL(i, :), angls);
    course = atan2(ned_gnss_vel(2), ned_gnss_vel(1)); % gnss angles

    % => Mechanisation
    g_tmp = comp_gravity(data.RawGNSS.LLA(i, :), "deg");
    if ~isnan(g_tmp)
        g = g_tmp;
    end
    omega = data.RawIMU.W(i, :) - bias_gyro;
    ac = a_centripetal(end, :).';
    vp = (Rbn * (data.RawIMU.SF(i, :).' * g - bias_acc.' - ac) + [0; 0; g]).';
    pp = mech_velocity(end, :);
    % pp(3) = -pp(3);
    Rbnp = Rbn * sq(omega);
    vb = Rbn.' * mech_velocity(end, :).';
    a_centripetal(end+1, :) = (sq(omega) * vb).';
    mech_position(end+1, :) = mech_position(end, :) + pp*dt;
    mech_velocity(end+1, :) = mech_velocity(end, :) + vp*dt;
    Rbn = norm_DCM(Rbn + Rbnp*dt);
    mech_angles(end+1, :) = dcm2ea(Rbn);
    
    % Kalman init
    Fvr = -sq(Rbn * (data.RawIMU.SF(i, :) * g - bias_acc).');
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
    Qk = 0.5 * dt * (Fi * G * Qc * G' + G * Qc * G' * Fi');

    % => Compute the Kalman gain
    K = P_next * H.' * inv(H * P_next * H.' + Rk);

    kss_array(:, :, i) = K;
    K_diag_hist(:, i) = sum(K(1:6, :)' .* eye(6))';

    % Calculation of pitch, roll and heading by integrating data from gyroscopes.
    % lecture 11, slide 18.
    om = data.RawIMU.W(i, :);
    g_phi = angles_gyro(end, 1);
    g_theta = angles_gyro(end, 2);
    g_R = [1,  sin(g_phi)*tan(g_theta),    cos(g_phi) * tan(g_theta);
        0,  cos(g_phi),               -sin(g_phi);
        0,  sin(g_phi)/cos(g_theta),    cos(g_phi)/cos(g_theta)];

    rpy_p = (g_R * om.').';
    angles_gyro(end+1, :) = angles_gyro(end, :) + rpy_p*dt;

    % -> zk compute
    if (isnan(course))
        % => Time update, predict
    	% -> Extrapolate the state
		xk(end+1, :) = x_next;
		Pk = P_next;

    	% -> Extrapolate the uncertainty
    	x_next = Fi * xk(end, :).';
    	P_next = Fi * Pk * Fi' + Qk;
		continue;
	end

	% GNSS angles
    if (norm(data.RawGNSS.VEL(i, :)) <= 0.1)
        angles_gnss(end+1, :) = nan;
	else
		angles_gnss(end+1, :) = course;
	end
	% If GNSS data was received:
    gnss_velocity(end+1, :) = ned_gnss_vel.';
    gps_track(end+1, :) = lla2ned(deg2rad(data.RawGNSS.LLA(i, :)).', lla0.');
    gps_track_lla(end+1, :) = data.RawGNSS.LLA(i, :).';

    % => Compute zk
    p_gnss = lla2ned(deg2rad(data.RawGNSS.LLA(i, :)).',lla0.');
    v_gnss = ned_gnss_vel;
    p_mech = mech_position(end, :).';
    v_mech = mech_velocity(end, :).';
    zk = [p_gnss - p_mech; v_gnss - v_mech];
    expected = H * x_next;
    inn_array(end + 1, :) = (zk - expected).';
    inn = (K * inn_array(end, :).').';
    % inn(1:6) = max(inn(1:6), -abs(zk - H * x_next)');
    % inn(1:6) = min(inn(1:6), abs(zk - H * x_next)');
	% 
    xk(end + 1, :) = (x_next + inn.').';
    xk_history(end + 1, :) = xk(end, :);
    bias_acc = bias_acc + xk(end, 10:12);
    biases_gyro = bias_gyro + xk(13:15);
    Rbn = (I + sq(xk(end, 7:9))) * Rbn;

    % P - symetric positive definite!;
    Pk = (I15 - K * H) * P_next * (I15 - K * H)' + K * Rk * K';
    
    % prediction
    xk(end, 7:end) = 0;
	x_next = Fi*xk(end, :).';
	P_next = Fi*Pk * Fi' + Qk;

end

kalman_position_lla = rad2deg(ned2lla(kalman_position, lla0));

% Plot 2D LLA position
figure(10);
geoplot(gps_track_lla(:, 1),gps_track_lla(:, 2), kalman_position_lla(:, 1), kalman_position_lla(:, 2))
grid on
geobasemap topographic
title("LLA 2d position (bad precision duew to convertion)")

% Plot 2D positionr
figure(1);
plot(gps_track(:, 2), gps_track(:, 1));
hold on
plot(kalman_position(:, 2), kalman_position(:, 1));
grid on
title("NED 2d position")
legend("GPS", "Kalman filter")

% Plot state vector history
figure(111);
for i = 1:15;
     subplot(5, 3, i);
     plot(xk_history(:, i));
     xlabel("Time, [s]");
     title("xk el: " + i);
     grid on
end


% Plot position in every axis
figure(11);
subplot(3, 1, 1)
plot((1:size(gps_track, 1)) * dt_gnss, gps_track(:, 1))
hold on
plot((1:size(kalman_position, 1)) * dt, kalman_position(:, 1))
grid on
title("Position")
legend("GPS", "Kalman filter")
subplot(3, 1, 2)
plot((1:size(gps_track, 1)) * dt_gnss, gps_track(:, 2))
hold on
plot((1:size(kalman_position, 1)) * dt, kalman_position(:, 2))
grid on
legend("GPS", "Kalman filter")
subplot(3, 1, 3)
plot((1:size(gps_track, 1)) * dt_gnss, gps_track(:, 3))
hold on
plot((1:size(kalman_position, 1)) * dt, kalman_position(:, 3))
grid on
legend("GPS", "Kalman filter")

% PLOT angles
figure(2)
subplot(3, 1, 1)
plot((1:size(angles_acc, 1)) * dt, angles_acc(:, 1))
hold on
plot((1:size(kalman_angles, 1)) * dt, kalman_angles(:, 1))
hold on
plot((1:size(angles_gyro, 1)) * dt, angles_gyro(:, 1))
hold on
plot((1:size(mech_angles, 1)) * dt, mech_angles(:, 1))
title("Angles")
xlabel("Time, [s]")
ylabel("angle [rad]")
grid on
legend("accelerometer roll", "kalman roll", "gyro roll")

subplot(3, 1, 2)
plot((1:size(angles_acc, 1)) * dt, angles_acc(:, 2))
hold on
plot((1:size(kalman_angles, 1)) * dt, kalman_angles(:, 2))
hold on
plot((1:size(angles_gyro, 1)) * dt, angles_gyro(:, 2))
hold on
plot((1:size(mech_angles, 1)) * dt, mech_angles(:, 2))
xlabel("Time, [s]")
ylabel("angle [rad]")
grid on
legend("accelerometer pitch", "kalman pitch", "gyro pitch")

subplot(3, 1, 3)
scatter((1:size(angles_gnss, 1)) * dt_gnss, angles_gnss, "*")
hold on
plot((1:size(kalman_angles, 1)) * dt, kalman_angles(:, 3))
hold on
plot((1:size(angles_gyro, 1)) * dt, angles_gyro(:, 3))
hold on
plot((1:size(mech_angles, 1)) * dt, mech_angles(:, 3))
xlabel("Time, [s]")
ylabel("angle [rad]")
grid on
legend("gnss course", "kalman yaw", "gyro yaw")

% PLOT velocities
figure(3)
subplot(3, 1, 1);
plot((1:size(gnss_velocity, 1)) * dt_gnss, gnss_velocity(:, 1))
hold on
plot((1:size(kalman_velocity, 1)) * dt, kalman_velocity(:, 1))
title("Velocity in NED")
ylabel("Velocity N, [m/s]")
xlabel("Time, [s]")
grid on
legend("gnss velocity", "kalman velocity")

subplot(3, 1, 2);
plot((1:size(gnss_velocity, 1)) * dt_gnss, gnss_velocity(:, 2))
hold on
plot((1:size(kalman_velocity, 1)) * dt, kalman_velocity(:, 2))
ylabel("Velocity E, [m/s]")
xlabel("Time, [s]")
grid on
legend("gnss velocity", "kalman velocity")

subplot(3, 1, 3);
plot((1:size(gnss_velocity, 1)) * dt_gnss, gnss_velocity(:, 3))
hold on
plot((1:size(kalman_velocity, 1)) * dt, kalman_velocity(:, 3))
ylabel("Velocity D, [m/s]")
xlabel("Time, [s]")
grid on
legend("gnss velocity", "kalman velocity")

% PLOT P matrix elements, diagonal
figure(4);
for i = 1:15;
     subplot(5, 3, i);
     plot(pss_array(:, i));
     xlabel("Time, [s]")
     title("P diag el: " + i)
     grid on
end

% PLOT innovation
figure(5);
for i = [1 3 5];
     subplot(3, 2, i);
     plot(inn_array(:, i));
     xlabel("Time, [s]")
     title("innovation (position): " + i)
     grid on
end
for i = [2 4 6];
     subplot(3, 2, i);
     plot(inn_array(:, i));
     xlabel("Time, [s]")
     title("innovation (velocity): " + i)
     grid on
end
