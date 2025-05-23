clear; clc;

% Datasets
datasets = {'L2Data4.csv', 'L2Data6.csv'};
titulos = {'L2Data4', 'L2Data6'};

% Parameters
dt = 0.01;
Q_kalman = 1e-3 * eye(6);
Q_ekf = 1e-4 * eye(5);  % agora com 5 estados
R_kalman = diag([7.6e-5, 8.22e-4, 108.90, 100.65]);
R_ekf = diag([7.6e-5, 8.22e-4]);

% Initialization
mae_phi_k = zeros(1,2);
mae_theta_k = zeros(1,2);
mae_phi_e = zeros(1,2);
mae_theta_e = zeros(1,2);

for d = 1:2
    dataset = datasets{d};
    titulo = titulos{d};

    data = readtable(dataset, 'VariableNamingRule','modify');
    acc_x = data.acc_x;
    acc_y = data.acc_y;
    acc_z = data.acc_z;

    p_m = data.gyro_x;
    q_m = data.gyro_y;
    r_m = data.gyro_z;

    phi_cf = deg2rad(data.stateEstimate_roll);
    theta_cf = deg2rad(data.stateEstimate_pitch);

    phi_m = atan2(acc_y, acc_z);
    theta_m = atan2(-acc_x, sqrt(acc_y.^2 + acc_z.^2));

    N = length(phi_m);
    t = (0:N-1) * dt;

    %% Linear KF 
    A = [0 0 1 0 0 0;
         0 0 0 1 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0];
    Ad = eye(6) + dt * A;
    C_k = [1 0 0 0 0 0;
           0 1 0 0 0 0;
           0 0 1 0 1 0;
           0 0 0 1 0 1];

    x_k = zeros(6, 1);
    P_k = 10 * eye(6);
    x_est_kalman = zeros(6, N);

    for k = 1:N
        y_k = [phi_m(k); theta_m(k); p_m(k); q_m(k)];
        x_pred = Ad * x_k;
        P_pred = Ad * P_k * Ad' + Q_kalman;
        K = P_pred * C_k' / (C_k * P_pred * C_k' + R_kalman);
        x_k = x_pred + K * (y_k - C_k * x_pred);
        P_k = (eye(6) - K * C_k) * P_pred;
        x_est_kalman(:, k) = x_k;
    end

    phi_est_kalman = rad2deg(x_est_kalman(1, :));
    theta_est_kalman = rad2deg(x_est_kalman(2, :));

    %% EKF
    x_e = zeros(5, 1);  % [phi, theta, bp, bq, br]
    P_e = 10 * eye(5);
    x_est_ekf = zeros(5, N);

    for k = 1:N
        phi = x_e(1); theta = x_e(2);
        bp = x_e(3); bq = x_e(4); br = x_e(5);

        % Corrected measurements
        p = p_m(k) - bp;
        q = q_m(k) - bq;
        r = r_m(k) - br;

        % Non linear model
        phi_dot = p + sin(phi)*tan(theta)*q + cos(phi)*tan(theta)*r;
        theta_dot = cos(phi)*q - sin(phi)*r;
        f = [phi_dot; theta_dot; 0; 0; 0];

        x_pred = x_e + f * dt;

        % Jacobian
        A = eye(5);
        A(1,1) = 1 + dt*(cos(phi)*tan(theta)*q - sin(phi)*tan(theta)*r);
        A(1,2) = dt*(sin(phi)/(cos(theta)^2)*q + cos(phi)/(cos(theta)^2)*r);
        A(1,3) = -dt;
        A(1,4) = -dt*sin(phi)*tan(theta);
        A(1,5) = -dt*cos(phi)*tan(theta);
        A(2,1) = -dt*(sin(phi)*q + cos(phi)*r);
        A(2,4) = -dt*cos(phi);
        A(2,5) = dt*sin(phi);

        P_pred = A * P_e * A' + Q_ekf;

        H = [1 0 0 0 0;
             0 1 0 0 0];
        y = [phi_m(k); theta_m(k)];
        y_hat = H * x_pred;
        S = H * P_pred * H' + R_ekf;
        K = P_pred * H' / S;
        x_e = x_pred + K * (y - y_hat);
        P_e = (eye(5) - K * H) * P_pred;
        x_est_ekf(:, k) = x_e;
    end

    phi_est_ekf = rad2deg(x_est_ekf(1, :));
    theta_est_ekf = rad2deg(x_est_ekf(2, :));
    phi_cf_deg = rad2deg(phi_cf);
    theta_cf_deg = rad2deg(theta_cf);

    % MAE's
    mae_phi_k(d) = mean(abs(phi_est_kalman - phi_cf_deg'));
    mae_theta_k(d) = mean(abs(theta_est_kalman - theta_cf_deg'));
    mae_phi_e(d) = mean(abs(phi_est_ekf - phi_cf_deg'));
    mae_theta_e(d) = mean(abs(theta_est_ekf - theta_cf_deg'));

    % Graficos
    figure('Name',['Comparacao EKF vs Kalman - ' titulo]);
    subplot(2,1,1);
    plot(t, phi_cf_deg, 'k', t, phi_est_kalman, 'r--', t, phi_est_ekf, 'b');
    legend('Crazyflie','Kalman Linear','EKF');
    ylabel('\phi (Roll) [°]');
    title(['Comparacao de Roll - ' titulo]);

    subplot(2,1,2);
    plot(t, theta_cf_deg, 'k', t, theta_est_kalman, 'r--', t, theta_est_ekf, 'b');
    legend('Crazyflie','Kalman Linear','EKF');
    ylabel('\theta (Pitch) [°]');
    xlabel('Tempo [s]');
    title(['Comparacao de Pitch - ' titulo]);
end

% Tabela de MAEs
fprintf('\n====== Tabela de MAE (Erro Médio Absoluto) ======\n');
fprintf('| Dataset   | Filtro         | MAE Roll (°) | MAE Pitch (°) |\n');
fprintf('|-----------|----------------|--------------|---------------|\n');
for i = 1:2
    fprintf('| %-9s | Kalman Linear  |    %6.2f     |     %6.2f     |\n', titulos{i}, mae_phi_k(i), mae_theta_k(i));
    fprintf('| %-9s | EKF            |    %6.2f     |     %6.2f     |\n', titulos{i}, mae_phi_e(i), mae_theta_e(i));
end
fprintf('=================================================\n');
