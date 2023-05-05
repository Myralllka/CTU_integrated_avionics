%% Load data from file
load("2022-05-13_CarData.mat");

initial_heading = 177; % deg
accel_data = RawIMU.SF;
gyro_data = RawIMU.W;
% gyro_data = gyro_data - mean(gyro_data(1:1000, :));
dt = 0.01;

%% Get angles from accelerometer
pitches_acc = zeros([1, length(accel_data)]);
rolls_acc = zeros([1, length(accel_data)]);
% Assmuing that the car is static at first and g is not the same everywhere
% on Earth
g = abs(mean(sum(accel_data(1:1000, :), 2)));

for i = 1:length(accel_data)
    longitudal = accel_data(i, 1);
    lateral = accel_data(i, 2);
    vertical = accel_data(i, 3);

    % Assuming that the vehicle is not accelerating and the only measured
    % force is the gravity
    pitches_acc(i) = asin(longitudal / g);
    rolls_acc(i) = -asin(lateral / g);

end

%% Integrate gyro data to get attitude in each time point
gyro_data_i = gyro_data - mean(gyro_data(1:1000, :));
dt = 0.01;
pitches_gyro = zeros([1, length(gyro_data)]);
rolls_gyro = zeros([1, length(gyro_data)]);
heading_gyro = zeros([1, length(gyro_data)]);
current_heading = initial_heading;

roll = mean(rolls_acc(1:1000));
pitch = mean(pitches_acc(1:1000));
yaw = initial_heading / 180 * pi;

for i = 1:length(gyro_data)
    w_x = gyro_data_i(i, 1);
    w_y = gyro_data_i(i, 2);
    w_z = gyro_data_i(i, 3);

    droll = w_x + w_y * sin(roll) * tan(pitch) + w_z * cos(roll) * tan(pitch);
    dpitch = w_y * cos(roll) - w_z * sin(roll);
    dyaw = w_y * sin(roll) / cos(pitch) + w_z * cos(roll) / cos(pitch);

    roll = roll + droll * dt;
    pitch = pitch + dpitch * dt;
    yaw = yaw + dyaw * dt;
    yaw = change_range_angle(yaw, 0, 0, 1);
    if yaw > pi
        yaw = yaw - 2 * pi;
    end
    
    rolls_gyro(i) = roll;
    pitches_gyro(i) = pitch;
    heading_gyro(i) = yaw;

end


%% Calculate and velocity of GPS in NED
pos_LLA = RawGNSS.LLA;
vel_ECEF = RawGNSS.VEL;
lla_ref = [pos_LLA(2, 1) / 180 * pi, pos_LLA(2, 2) / 180 * pi, pos_LLA(2, 3)];

gps_pos_ned = [];
gps_vel_ned = [];
gps_pos_ned_filtered = [];
gps_vel_ned_filtered = [];

for i = 1:length(pos_LLA)
    if (isnan(pos_LLA(i, 1)))
        gps_pos_ned(:, end + 1) = [NaN; NaN; NaN];
        gps_vel_ned(:, end + 1) = [NaN; NaN; NaN];
        continue;
    end
    lat = pos_LLA(i, 1);
    lon = pos_LLA(i, 2);
    a = pos_LLA(i, 3);
    gps_pos_ned(:, end + 1) = lla2ned([lat / 180 * pi, lon / 180 * pi, a]', lla_ref');
    gps_vel_ned(:, end + 1) = ecef2ned(RawGNSS.LLA(i, :), vel_ECEF(i, :));
    gps_pos_ned_filtered(:, end + 1) = gps_pos_ned(:, end);
    gps_vel_ned_filtered(:, end + 1) = gps_vel_ned(:, end);

end
   


%% Kalman filter part

% Initialization of mechanization variables
initial_roll = mean(rolls_acc(1:1000));
initial_pitch = mean(pitches_acc(1:1000));
g = comp_gravity(RawGNSS.LLA(2, :), "deg");
g_n = [0; 0; g];
v_mech = [0; 0; 0];
pos_mech = [0; 0; 0];
a_c = [0; 0; 0];

LLA_rad = RawGNSS.LLA;
LLA_rad(:, 1) = deg2rad(LLA_rad(:, 1));
LLA_rad(:, 2) = deg2rad(LLA_rad(:, 2));

ref_LLA = LLA_rad(2, :)';

% Initialization of Kalman filter variables
H = [eye(6), zeros([6, 9])];

% Define Fs
F_pp = zeros(3);
F_pro = zeros(3);
F_pv = eye(3);
F_pv(3, 3) = 1;

F_vp = zeros(3);
F_vv = zeros(3);
F_va = eye(3);

F_rop = zeros(3);
F_rov = zeros(3);
F_roro = zeros(3);
F_rog = eye(3);

% Random constant bias 
F_aa = zeros(3);
F_gg = zeros(3);


x_init = zeros([15, 0]);
x_init(7, 1) = initial_roll;
x_init(8, 1) = initial_pitch;
x_init(9, 1) = initial_heading / 180 * pi;

x_init(10:12) = 0;
x_init(13:15) = mean(RawIMU.W(1:1000, :));


x = x_init;
x_p = x_init;
dcm = ea2dcm([initial_roll, initial_pitch, initial_heading / 180 * pi]);
dcm_hat = dcm;

% Find out the acceleration multipllier considering that in the beginning
% it is static and should measure exactly 1g
acc_mult = abs(mean(sum(accel_data(1:1000, :), 2)));


% TODO: find some good values for this
Q = eye(15);
Q(1:3, 1:3) = eye(3) * 0.1;
Q(4:6, 4:6) = diag(var(accel_data(1:1000, :))) * 10000;
Q(7:9, 7:9) = diag(var(gyro_data(1:1000, :))) * 10000;
Q(10:15, 10:15) = eye(6) * 0.0001;

R_c = eye(6);

% R_c(1:3, 1:3) = diag(var(gps_pos_ned_filtered(:, 1:100)'));
% R_c(4:6, 4:6) = diag(var(gps_vel_ned_filtered(:, 1:100)'));
R_c(1:3, 1:3) = eye(3) * 1.5^2;
R_c(4:6, 4:6) = eye(3) * 0.01^2;
R_k = R_c;


% Generate an initial matrix A
A = [F_pp, F_pv, F_pro, zeros([3, 6]);
     F_vp, F_vv, zeros([3, 9]);
     F_rop, F_rov, F_roro, zeros([3, 6]);
     zeros([3, 9]), F_aa, zeros(3);
     zeros([3, 12]), F_gg];

state_vectors = [x_init];

P_k = eye(15) * 50;
P_p = P_k;

N = length(accel_data);

% Initializa matrices for saving values history
pos_estimate = zeros([3, N]);
vel_estimate = zeros([3, N]);
P_hist = zeros([15, N]);
gyro_bias_hist = zeros([3, N]);
accel_bias_hist = zeros([3, N]);
accel_bias = [0; 0; 0];
gyro_bias = [0; 0; 0];
K_diag_hist = zeros([6, N]);
angles_hist = zeros([3, N]);


for i = 1:N
    if (isnan(x(1)))
        disp(i);
        break;
    end

    % Save values history to corresponding arrays
    pos_estimate(:, i) = pos_mech + x(1:3);
    vel_estimate(:, i) = v_mech + x(4:6);
    P_hist(:, i) = sum(P_k .* eye(15))';
    gyro_bias_hist(:, i) = gyro_bias;
    accel_bias_hist(:, i) = accel_bias;
    a = dcm2ea(dcm);
    phi = a(1);
    theta = a(2);
    psi = a(3);
    angles_hist(:, i) = [phi; theta; psi];


    % Retrieve sensor values, transform force to m/s^2 and subtract biases
    f_b = accel_data(i, :)' / acc_mult * g - accel_bias;
    w_b = gyro_data(i, :)' - gyro_bias;


    % Perform mechanization
    d_v = dcm * (f_b - a_c) + g_n;
    d_dcm = dcm * sq(w_b);
    v_mech = v_mech + dt * d_v;
    pos_mech = pos_mech + v_mech * dt;
    v_b = dcm \ v_mech;
    a_c = cross(w_b, v_b);
    dcm = norm_DCM(dcm + dt * d_dcm);
    % END MECHANIZATION


    % Update matrix A for current step
    f_n = dcm * f_b;

    F_vro = [0, f_n(3), -f_n(2); 
            -f_n(3), 0, f_n(1);
            f_n(2), -f_n(1), 0];
    A(4:6, 7:9) = -F_vro;
    A(4:6, 10:12) = -dcm;
    A(7:9, 13:15) = dcm;
    
    % Calculate model discretizatoin
    PHI = eye(15) + A * dt; 

    G = [eye(3), zeros([3, 12]);
     zeros(3), -dcm, zeros([3, 9]);
     zeros([3, 6]), dcm, zeros([3, 6]);
     zeros([6, 9]), eye(6)];
    Q_k = 0.5 * dt * (PHI * G * Q * G' + G * Q * G' * PHI');
    K_k = P_p * H' * inv(H * P_p * H' + R_k);

    % Save K
    K_diag_hist(:, i) = sum(K_k(1:6, :)' .* eye(6))';


    % If no GPS measurement, set values of x and P to the predicted ones
    % and predict new values
    if (isnan(RawGNSS.LLA(i, 1)))
        x = x_p;
        P_k = P_p;

        x_p = PHI * x;
        P_p = PHI * P_k * PHI' + Q_k;

        continue;
    end

    % When GNSS data is available, perform correctoin step
    pos_LLA = LLA_rad(i, :)';
    gps_pos = lla2ned(pos_LLA, ref_LLA);
    gps_vel = ecef2ned(RawGNSS.LLA(i, :), RawGNSS.VEL(i, :));
    
    z = [gps_pos; gps_vel] - [pos_mech; v_mech];

    inn = K_k * (z - H * x_p);

    % Scale innovation for p and v so its never greater than 1
    inn(1:6) = max(inn(1:6), -abs(z - H * x_p));
    inn(1:6) = min(inn(1:6), abs(z - H * x_p));

    % Update state vector and biases
    x = x_p + inn;
    gyro_bias = gyro_bias + x(13:15);
    accel_bias = accel_bias + x(10:12);
    dcm = norm_DCM(dcm - sq(x(7:9)));
    x(7:end) = 0;

    % Predict new x and P
    P_k = (eye(15) - K_k * H) * P_p * (eye(15) - K_k * H)' + K_k * R_k * K_k';
    x_p = PHI * x;
    P_p = PHI * P_k * PHI' + Q_k; 
end



%% Plot position in 2d


plot(pos_estimate(1, :), pos_estimate(2, :));
hold on;
plot(gps_pos_ned_filtered(1, :), gps_pos_ned_filtered(2, :));
hold off;
title("2d Position");
legend(["KF estimation", "GNSS"]);

%% Plot each position component on different subplots
figure
subplot(3, 1, 1);
plot(pos_estimate(1, :));
hold on;
plot(1:10:(10 * length(gps_pos_ned_filtered)), gps_pos_ned_filtered(1, :));
title("Position in NED (N)");
legend(["Estimation", "GNSS"]);

subplot(3, 1, 2);
plot(pos_estimate(2, :));
hold on;
plot(1:10:(10 * length(gps_pos_ned_filtered)), gps_pos_ned_filtered(2, :));
title("Position in NED (E)");
legend(["Estimation", "GNSS"]);


subplot(3, 1, 3);
plot(pos_estimate(3, :));
hold on;
plot(1:10:(10 * length(gps_pos_ned_filtered)), gps_pos_ned_filtered(3, :));
title("Position in NED (D)");
legend(["Estimation", "GNSS"]);


%% Plot velocity

figure
subplot(3, 1, 1);
plot(vel_estimate(1, :));
hold on;
plot(1:10:(10 * length(gps_vel_ned_filtered)), gps_vel_ned_filtered(1, :));
title("Velocity in NED (N)");
legend(["Estimation", "GNSS"]);

subplot(3, 1, 2);
plot(vel_estimate(2, :));
hold on;
plot(1:10:(10 * length(gps_vel_ned_filtered)), gps_vel_ned_filtered(2, :));
title("Velocity in NED (E)");
legend(["Estimation", "GNSS"]);


subplot(3, 1, 3);
plot(vel_estimate(3, :));
hold on;
plot(1:10:(10 * length(gps_vel_ned_filtered)), gps_vel_ned_filtered(3, :));
title("Velocity in NED (D)");
legend(["Estimation", "GNSS"]);


%% Plot angles

figure
subplot(3, 1, 1);
plot(rolls_acc);
hold on
plot(angles_hist(1, :)');
plot(rolls_gyro);
legend("Roll accelerometer", "Roll KF", "Roll gyro")

subplot(3, 1, 2);
plot(pitches_acc);
hold on
plot(angles_hist(2, :)');
plot(pitches_gyro);
legend("Pitch accelerometer", "Pitch KF", "Pitch gyro")

subplot(3, 1, 3);
hold on
plot(angles_hist(3, :)');
plot(heading_gyro);
legend("Yaw KF", "Yaw gyro");


%%
%% Matrix P
figure
subplot(3, 3, 1);
plot(P_hist(1, :));
title("P[1]");

subplot(3, 3, 2);
plot(P_hist(2, :));
title("P[2]");

subplot(3, 3, 3);
plot(P_hist(3, :));
title("P[3]");

subplot(3, 3, 4);
plot(P_hist(4, :));
title("P[4]");

subplot(3, 3, 5);
plot(P_hist(5, :));
title("P[5]");

subplot(3, 3, 6);
plot(P_hist(6, :));
title("P[6]");

subplot(3, 3, 7);
plot(P_hist(7, :));
title("P[7]");

subplot(3, 3, 8);
plot(P_hist(8, :));
title("P[8]");

subplot(3, 3, 9);
plot(P_hist(9, :));
title("P[9]");