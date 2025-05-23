clear; clc;

% Leitura dos datasets
datasets = {'L2Data4.csv', 'L2Data6.csv'};
titulos = {'L2Data4', 'L2Data6'};

% Definição dos parâmetros
dt = 0.01;
Q = 1e-3 * eye(6);
R = diag([7.6e-5, 4.7e-5, 108.90, 100.65]);

A = [0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];
Ad = eye(6) + dt * A;

C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 1 0;
     0 0 0 1 0 1];

% Inicialização
mae_phi_all = zeros(1, 2);
mae_theta_all = zeros(1, 2);

for i = 1:length(datasets)
    % Processamento dos dados
    filename = datasets{i};
    data = readtable(filename, 'VariableNamingRule', 'modify');

    acc_x = data.acc_x;
    acc_y = data.acc_y;
    acc_z = data.acc_z;

    gyro_x = data.gyro_x;
    gyro_y = data.gyro_y;

    % Graus para radianos
    phi_cf = deg2rad(data.stateEstimate_roll);
    theta_cf = deg2rad(data.stateEstimate_pitch);

    % Estimativa dos ângulos
    phi_m = atan2(acc_y, acc_z);
    theta_m = atan2(-acc_x, sqrt(acc_y.^2 + acc_z.^2));

    p_m = gyro_x;
    q_m = gyro_y;

    % Debug 
    fprintf('\n[%s] Verificação de variância:\n', titulos{i});
    fprintf('Var(phi_m): %.6f | Var(phi_cf): %.6f\n', var(phi_m), var(phi_cf));
    fprintf('Var(p_m): %.6f | Var(gyro_x): %.6f\n', var(p_m), var(gyro_x));

    N = length(phi_m);
    x = zeros(6, 1);
    P = 10 * eye(6);
    x_est = zeros(6, N);

    % KF
    for k = 1:N
        y = [phi_m(k); theta_m(k); p_m(k); q_m(k)];
        x_pred = Ad * x;
        P_pred = Ad * P * Ad' + Q;
        K = P_pred * C' / (C * P_pred * C' + R);
        x = x_pred + K * (y - C * x_pred);
        P = (eye(6) - K * C) * P_pred;
        x_est(:, k) = x;

        % Debug a cada 500 amostras
        if mod(k, 500) == 0
            fprintf('Amostra %d: phi estimado = %.4f rad\n', k, x(1));
        end
    end

    % Resultados
    phi_est = x_est(1, :);
    theta_est = x_est(2, :);

    % Converção para graus 
    phi_est_deg = rad2deg(phi_est);
    theta_est_deg = rad2deg(theta_est);
    phi_cf_deg = rad2deg(phi_cf);
    theta_cf_deg = rad2deg(theta_cf);

    mae_phi = mean(abs(phi_est_deg - phi_cf_deg'));
    mae_theta = mean(abs(theta_est_deg - theta_cf_deg'));

    mae_phi_all(i) = mae_phi;
    mae_theta_all(i) = mae_theta;

    fprintf('[%s] MAE Roll: %.2f° | MAE Pitch: %.2f°\n', titulos{i}, mae_phi, mae_theta);

    % Gráficos
    figure('Name', titulos{i});
    subplot(2,1,1);
    plot(phi_cf_deg, 'k', 'DisplayName','Crazyflie');
    hold on;
    plot(phi_est_deg, 'r', 'DisplayName','Kalman');
    ylabel('\phi (Roll) [°]');
    legend;
    title([titulos{i} ' - Estimativa de Roll']);

    subplot(2,1,2);
    plot(theta_cf_deg, 'k', 'DisplayName','Crazyflie');
    hold on;
    plot(theta_est_deg, 'b', 'DisplayName','Kalman');
    ylabel('\theta (Pitch) [°]');
    xlabel('Amostras');
    legend;
    title([titulos{i} ' - Estimativa de Pitch']);
end

% Resultados
fprintf('\n====== Tabela de MAE (Erro Médio Absoluto) ======\n');
fprintf('| Dataset   | MAE Roll (°) | MAE Pitch (°) |\n');
fprintf('|-----------|--------------|---------------|\n');
for i = 1:length(datasets)
    fprintf('| %-9s |    %6.2f     |     %6.2f     |\n', titulos{i}, mae_phi_all(i), mae_theta_all(i));
end
fprintf('===============================================\n');
