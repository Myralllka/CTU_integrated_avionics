clear all;
close all;

% Initial alignment, Mechanization of navigation equations

% 1) Implement a navigation equation mechanization algorithm that 
% will process data from file. Implement functions for:
% a. done - Conversion between the ECEF and NED systems (ecef2ned),
% b. done - Conversion of Euler angles to DCM (ea2dcm),
% c. done - Conversion of DCM to Euler angles (dcm2ea),
% d. given - Function for normalization of the DCM matrix (normDCM)
% 2) Next, implement the calculation of position angles according to the following points:
% a. done - Calculation of pitch and roll from triaxial accelerometer data,
% b. done - Course calculation from GNSS data,
% c. done - Calculation of pitch, roll and heading by integrating data from gyroscopes.
% 3) Implement the conversion of GNSS data from the Latitude-Longitude-Altitude (LLA) system to the system
% North-East-Down (NED).
% 4) Plot the following variables:
% a. Position angles obtained from mechanization, from accelerometer, gyroscope, GNSS,
% b. Position calculated by mechanization in NED,
% c. Speed calculated in NED.

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

roll_init_acc = mean(atan(data.RawIMU.SF(1:1000, 2) ./ data.RawIMU.SF(1:1000, 3))); % acc angles
pitch_init_acc = mean(-atan(data.RawIMU.SF(1:1000, 1) ./ data.RawIMU.SF(1:1000, 3)));

init_angle = [roll_init_acc, pitch_init_acc, deg2rad(177)];
C = ea2dcm(init_angle);

position = [[0, 0, 0]];
velocity = [[0, 0, 0]];
a_centripetal = [[0, 0, 0]];

angles_acc = [[roll_init_acc, pitch_init_acc]];
angles_gyro = [[roll_init_acc, pitch_init_acc, init_angle(3)]];
angles_gnss = [[init_angle(3)]];

angles_mechanisation = [[0, 0, 0]];

lla0 = deg2rad(data.RawGNSS.LLA(2, :));
gps_track = [];
N = size(data.RawGNSS.LLA, 1);

for i = 1001:N;

    % Calculation of pitch and roll from triaxial accelerometer data
    roll = atan(data.RawIMU.SF(i, 2) / data.RawIMU.SF(i, 3)); % acc angles
    pitch = -atan(data.RawIMU.SF(i, 1) / data.RawIMU.SF(i, 3));
    angles_acc(end+1, :) = [roll, pitch];

    % Course calculation from GNSS data
    angls = [deg2rad(data.RawGNSS.LLA(i, 1)), deg2rad(data.RawGNSS.LLA(i, 2)), data.RawGNSS.LLA(i, 3)];
    ned_gnss_vel = ecef2ned(data.RawGNSS.VEL(i, :), angls);
    course = atan2(ned_gnss_vel(2), ned_gnss_vel(1)); % gnss angles
    if ~isnan(course)
        gps_track(end+1, :) = lla2ned(deg2rad(data.RawGNSS.LLA(i, :)).', lla0.');
        if (norm(data.RawGNSS.VEL(i, :)) <= 0.1)
            angles_gnss(end+1, :) = nan;
        else
            angles_gnss(end+1, :) = course;
        end
    end

    % Mechanisation
    g_tmp = comp_gravity(data.RawGNSS.LLA(i, :), "deg");
    if ~isnan(g_tmp)
        g = g_tmp;
    end

    omega = data.RawIMU.W(i, :);
    ac = a_centripetal(end, :).';
    vp = (C * (data.RawIMU.SF(i, :).' * g - ac) + [0; 0; g]).';
    pp = velocity(end, :);
    pp(3) = -pp(3);
    Cp = C * sq(omega);

    vb = C.' * velocity(end, :).';
    a_centripetal(end+1, :) = (sq(omega) * vb).';

    position(end+1, :) = position(end, :) + pp*dt;
    velocity(end+1, :) = velocity(end, :) + vp*dt;
    C = norm_DCM(C + Cp*dt);

    angles_mechanisation(end+1, :) = dcm2ea(C);

    % Calculation of pitch, roll and heading by integrating data from gyroscopes.
    % lecture 11, slide 18.
    phi = angles_gyro(end, 1);
    theta = angles_gyro(end, 2);
    R = [1,  sin(phi)*tan(theta),    cos(phi) * tan(theta);
         0,  cos(phi),               -sin(phi);
         0,  sin(phi)/cos(theta),    cos(phi)/cos(theta)];

    om = data.RawIMU.W(i, :);
    rpy_p = (R * om.').';
    angles_gyro(end+1, :) = angles_gyro(end, :) + rpy_p*dt;

end

figure(1);

subplot(3, 1, 1)
plot((1:size(angles_acc, 1)) * dt, angles_acc(:, 1))
hold on
plot((1:size(angles_mechanisation, 1)) * dt, angles_mechanisation(:, 1))
hold on
plot((1:size(angles_gyro, 1)) * dt, angles_gyro(:, 1))
title("Angles")
xlabel("Time, [s]")
ylabel("angle [rad]")
grid on
legend("accelerometer roll", "mechanisation roll", "gyro roll")

subplot(3, 1, 2)
plot((1:size(angles_acc, 1)) * dt, angles_acc(:, 2))
hold on
plot((1:size(angles_mechanisation, 1)) * dt, angles_mechanisation(:, 2))
hold on
plot((1:size(angles_gyro, 1)) * dt, angles_gyro(:, 2))
xlabel("Time, [s]")
ylabel("angle [rad]")
grid on
legend("accelerometer pitch", "mechanisation pitch", "gyro pitch")

subplot(3, 1, 3)
scatter((1:size(angles_gnss, 1)) * dt_gnss, angles_gnss, "*")
hold on
plot((1:size(angles_mechanisation, 1)) * dt, angles_mechanisation(:, 3))
hold on
plot((1:size(angles_gyro, 1)) * dt, angles_gyro(:, 3))
xlabel("Time, [s]")
ylabel("angle [rad]")
grid on
legend("gnss course", "mechanisation yaw", "gyro yaw")


figure(2);
subplot(3, 1, 1);
plot(velocity(:, 1));
title("Velocity in NED")
ylabel("Velocity N, [m/s]")
xlabel("Time, [s]")
grid on

subplot(3, 1, 2);
plot(velocity(:, 2));
ylabel("Velocity E, [m/s]")
xlabel("Time, [s]")
grid on

subplot(3, 1, 3);
plot(velocity(:, 3));
ylabel("Velocity D, [m/s]")
xlabel("Time, [s]")
grid on

figure(3);
subplot(3, 1, 1);
plot(position(:, 1));
xlabel("Time, [s]")
ylabel("Position N, [m]")
title("Position, NED")
grid on

subplot(3, 1, 2);
plot(position(:, 2));
xlabel("Time, [s]")
ylabel("Position N, [m]")
grid on

subplot(3, 1, 3);
plot(position(:, 3));
xlabel("Time, [s]")
ylabel("Position N, [m]")
grid on

figure(4)
plot(gps_track(:, 2), gps_track(:, 1));
hold on
plot(position(:, 2), position(:, 1));
grid on