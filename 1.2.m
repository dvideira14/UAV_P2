data2 = readtable('L2Data2.csv');
n = height(data2);
hover_data = data2(round(n*0.35):round(n*0.65), :);
figure;
subplot(2,1,1);
plot(hover_data.timestamp, hover_data.acc_x, hover_data.timestamp, hover_data.acc_y, hover_data.timestamp, hover_data.acc_z);
legend('acc.x','acc.y','acc.z');
title('Acelerómetro - Hover (L2Data2)');
xlabel('Tempo [ms]');
ylabel('Aceleração [m/s^2]');
subplot(2,1,2);
plot(hover_data.timestamp, hover_data.gyro_x, hover_data.timestamp, hover_data.gyro_y, hover_data.timestamp, hover_data.gyro_z);
legend('gyro.x','gyro.y','gyro.z');
title('Giroscópio - Hover (L2Data2)');
xlabel('Tempo [ms]');
ylabel('Velocidade angular [rad/s]');

saveas(gcf, 'L2Data2_Hover_Sensores.png');

acc_hover_mean = mean([hover_data.acc_x, hover_data.acc_y, hover_data.acc_z]);
acc_hover_std = std([hover_data.acc_x, hover_data.acc_y, hover_data.acc_z]);
acc_hover_var = var([hover_data.acc_x, hover_data.acc_y, hover_data.acc_z]);

gyro_hover_mean = mean([hover_data.gyro_x, hover_data.gyro_y, hover_data.gyro_z]);
gyro_hover_std = std([hover_data.gyro_x, hover_data.gyro_y, hover_data.gyro_z]);
gyro_hover_var = var([hover_data.gyro_x, hover_data.gyro_y, hover_data.gyro_z]);

disp('Acelerómetro (Hover - L2Data2)');
disp(acc_hover_mean); disp(acc_hover_std); disp(acc_hover_var);

disp('Giroscópio (Hover - L2Data2)');
disp(gyro_hover_mean); disp(gyro_hover_std); disp(gyro_hover_var);
