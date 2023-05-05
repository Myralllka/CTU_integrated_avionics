clear all;
close all;
% Create an matlab-file to demonstrate Kalman Filtering - 
% Constant estimation. Matlab-file will have to meet the following points:
% done - a) Generate a random signal with a defined standard deviation of length 2000
% samples, the signal will contain a step change in level, i.e. first 1000 level samples
% L1, second 1000 samples of level L2.
% done - b) Implement the KF to estimate the constant and tune the filter properly 
% according to the instructions given at the lecture.
% done - c) Plot the KF parameters:
% c1) Raw signal unfiltered vs. state vector (x);
% c2) Innovation;
% c3) Kalman gain;
% c4) Error covariance matrix P.
constant1 = 10;
constant2 = 15;
N = 2000;
dt = 0.01;

x = zeros(1, N);
x(1:1000) = x(1:1000) + constant1;
x(1001:end) = x(1001:end) + constant2;
noise = randn(1, N);
sig = x + noise;

% Values initialisation
A = [0];

xk = [0];

Pk = [1];
Qk = [std(sig)^2];

R = std(sig)^2 / dt;

G = [1];
I = eye(1);
Fi = [1];
H = [1];
inn_array = [0];
pss_array = [Pk];
kss_array = [0];

% x-, Pk- - extrapolation
% x+, Pk+ - observation

for i = 1:N
    % time step
    x_next_extra = Fi * xk(i);
    P_next_extra = Fi * Pk * Fi.' + Qk;
    
    % 
    K = P_next_extra * H.' / (H * P_next_extra * H.'  + R);
    
    zk = sig(i);

    expected = H * x_next_extra;
    inn_array(i) = zk - expected;

    xk(i + 1) = x_next_extra + K * inn_array(i);

    % P - symetric positive definite!
    Pk = (I - K*H) * P_next_extra;
    Pk = (Pk + Pk.') / 2;
    pss_array(i) = Pk;
    kss_array(i) = K;
    % Pk = (I - K*H)*P_next * (I - K*H).' + K*R*K.';

    %
    % diag(Pk) - dispersion of values that are looked
end


figure(1);

subplot(4, 1, 1);
grid on
plot(sig, '-g')
hold on
plot(xk(1:end-1), '-b')
title('Original signal (green) and filtered signal(blue).');

%
subplot(4, 1, 2);
plot(inn_array, '-r')
grid on
title('Innovation.');

%
subplot(4, 1, 3);
plot(kss_array, '-r')
grid on
title('Kalman gain.');

% 
subplot(4, 1, 4);
plot(pss_array, '-r')
grid on
title('Error covariance matrix P.');
