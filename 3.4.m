% Comparação EKF com e sem Bias (L2Data6)
clear; clc;

dataset = 'L2Data6.csv';
dt = 0.01;
g = 9.81;

% Covariâncias
Q = 1e-4 * eye(5);     % com bias
Q2 = 1e-4 * eye(2);    % sem bias
R = diag([7.6e-5, 4.7e-5]);

data = readtable(dataset, 'VariableNamingRule','modify');

% Dados dos sensores
acc_x = data.acc_x;
acc_y = data.acc_y;
acc_z = data.acc_z;
p_m = data.gyro_x;
q_m = data.gyro_y;
phi_cf = deg2rad(data.stateEstimate_roll);
theta_cf = deg2rad(data.stateEstimate_pitch);

% Inclinômetro
phi_m = atan2(acc_y, acc_z);
theta_m = atan2(-acc_x, sqrt(acc_y.^2 + acc_z.^2));
z = [phi_m'; theta_m'];

N = length(phi_m);
t = (0:N-1) * dt;

%% EKF com bias (5 estados: phi, theta, bp, bq, br)
x_e = zeros(5,1); P = eye(5);
x_est_bias = zeros(5,N);

for k = 2:N
    phi = x_e(1); theta = x_e(2);
    bp = x_e(3); bq = x_e(4);

    p = p_m(k) - bp;
    q = q_m(k) - bq;

    % Dinâmica não linear
    phi_dot = p + sin(phi)*tan(theta)*q;
    theta_dot = cos(phi)*q;
    f = [phi_dot; theta_dot; 0; 0; 0];
    x_pred = x_e + f * dt;

    % Jacobiano
    A = eye(5);
    A(1,1) = 1 + dt * cos(phi)*tan(theta)*q;
    A(1,2) = dt * (sin(phi)/(cos(theta)^2)) * q;
    A(1,3) = -dt;
    A(1,4) = -dt * sin(phi)*tan(theta);
    A(2,1) = -dt * sin(phi)*q;
    A(2,4) = -dt * cos(phi);

    P_pred = A * P * A' + Q;

    H = [1 0 0 0 0;
         0 1 0 0 0];
    y = z(:,k);
    y_hat = H * x_pred;
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;

    x_e = x_pred + K * (y - y_hat);
    P = (eye(5) - K * H) * P_pred;
    x_est_bias(:,k) = x_e;
end

%% EKF sem bias (2 estados: phi, theta)
x_e2 = zeros(2,1); P2 = eye(2);
x_est_nobias = zeros(2,N);

for k = 2:N
    phi = x_e2(1); theta = x_e2(2);
    p = p_m(k); q = q_m(k);

    phi_dot = p + sin(phi)*tan(theta)*q;
    theta_dot = cos(phi)*q;
    f = [phi_dot; theta_dot];
    x_pred = x_e2 + f * dt;

    A2 = eye(2);
    A2(1,1) = 1 + dt * cos(phi)*tan(theta)*q;
    A2(1,2) = dt * (sin(phi)/(cos(theta)^2)) * q;
    A2(2,1) = -dt * sin(phi) * q;

    P2_pred = A2 * P2 * A2' + Q2;
    H2 = eye(2);
    y = z(:,k);
    y_hat = H2 * x_pred;
    S2 = H2 * P2_pred * H2' + R;
    K2 = P2_pred * H2' / S2;

    x_e2 = x_pred + K2 * (y - y_hat);
    P2 = (eye(2) - K2 * H2) * P2_pred;
    x_est_nobias(:,k) = x_e2;
end

% Converter para graus
phi_bias = rad2deg(x_est_bias(1,:));
theta_bias = rad2deg(x_est_bias(2,:));
phi_nobias = rad2deg(x_est_nobias(1,:));
theta_nobias = rad2deg(x_est_nobias(2,:));
phi_cf = rad2deg(phi_cf);
theta_cf = rad2deg(theta_cf);

%% Gráfico comparativo
figure('Name','EKF com e sem Bias - L2Data6','NumberTitle','off');
subplot(2,1,1);
plot(t, phi_cf, 'k--', t, phi_bias, 'b', t, phi_nobias, 'r');
xlabel('Tempo [s]'); ylabel('\\phi (Roll) [°]');
title('Comparação de Roll - EKF com e sem Bias');
legend('Crazyflie','EKF com bias','EKF sem bias'); grid on;

subplot(2,1,2);
plot(t, theta_cf, 'k--', t, theta_bias, 'b', t, theta_nobias, 'r');
xlabel('Tempo [s]'); ylabel('\\theta (Pitch) [°]');
title('Comparação de Pitch - EKF com e sem Bias');
legend('Crazyflie','EKF com bias','EKF sem bias'); grid on;
