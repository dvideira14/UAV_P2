% Carregar dados do datset
data = readtable('L2Data4.csv');

% Parâmetros
dt = 0.01; % passo de amostragem [s]
N = height(data);

% Extract measurements
acc = [data.acc_x, data.acc_y, data.acc_z];
gyro = [data.gyro_x, data.gyro_y]; % [p, q]
roll_true = deg2rad(data.stateEstimate_roll);
pitch_true = deg2rad(data.stateEstimate_pitch);

% Inicializações
x_with_bias = zeros(N, 4);    % [phi, theta, bp, bq]
x_no_bias = zeros(N, 2);      % [phi, theta]

Pwb = eye(4) * 0.1;           % Covariância inicial (com bias)
Pnb = eye(2) * 0.1;           % Covariância inicial (sem bias)

Qwb = eye(4) * 0.001;         % Ruído de processo ajustado (com bias)
Qnb = eye(2) * 0.001;         % Ruído de processo ajustado (sem bias)

R = eye(2) * 0.05;            % Ruído de medição (ajustado)

% Loop principal
for k = 2:N
    % Medições atuais
    p = gyro(k,1);
    q = gyro(k,2);
    ax = acc(k,1);
    ay = acc(k,2);
    az = acc(k,3);

    % Estimativas dos ângulos a partir do acelerómetro
    phi_acc = atan2(ay, az);
    theta_acc = atan2(-ax, sqrt(ay^2 + az^2));
    z = [phi_acc; theta_acc];

    %% With bias
    x_prev = x_with_bias(k-1,:)';

    % Previsão do estado
    phi_pred   = x_prev(1) + dt * (p - x_prev(3));
    theta_pred = x_prev(2) + dt * (q - x_prev(4));
    x_pred = [phi_pred; theta_pred; x_prev(3); x_prev(4)];

    % Matriz A
    A = [1, 0, -dt,  0;
         0, 1,  0, -dt;
         0, 0,  1,   0;
         0, 0,  0,   1];

    % Previsão da covariância
    P_pred = A * Pwb * A' + Qwb;

    % Medição esperada
    H = [1, 0, 0, 0;
         0, 1, 0, 0];

    % Inovação e ganho
    y = z - H * x_pred;
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;

    % Actualização
    x_with_bias(k,:) = (x_pred + K * y)';
    Pwb = (eye(4) - K * H) * P_pred;

    %% without bias
    x_prev_nb = x_no_bias(k-1,:)';

    % Previsão do estado
    x_pred_nb = x_prev_nb + dt * [p; q];

    % Previsão da covariância
    P_pred_nb = Pnb + Qnb;

    % Medição esperada
    H_nb = eye(2);

    % Gain
    y_nb = z - H_nb * x_pred_nb;
    S_nb = H_nb * P_pred_nb * H_nb' + R;
    K_nb = P_pred_nb * H_nb' / S_nb;

    % Actualização
    x_no_bias(k,:) = (x_pred_nb + K_nb * y_nb)';
    Pnb = (eye(2) - K_nb * H_nb) * P_pred_nb;
end

%% Comparação
tempo = (0:N-1) * dt;

figure;
subplot(2,1,1);
plot(tempo, roll_true, 'k', 'LineWidth', 1.2); hold on;
plot(tempo, x_with_bias(:,1), 'b--', 'LineWidth', 1.2);
plot(tempo, x_no_bias(:,1), 'r:', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Roll [rad]');
legend('Real', 'Com bias', 'Sem bias');
title('Estimativa do ângulo de Roll');

subplot(2,1,2);
plot(tempo, pitch_true, 'k', 'LineWidth', 1.2); hold on;
plot(tempo, x_with_bias(:,2), 'b--', 'LineWidth', 1.2);
plot(tempo, x_no_bias(:,2), 'r:', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Pitch [rad]');
legend('Real', 'Com bias', 'Sem bias');
title('Estimativa do ângulo de Pitch');
