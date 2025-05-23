data1 = readtable('L2Data1.csv');

figure;
subplot(2,1,1);
plot(data1.timestamp, data1.acc_x, data1.timestamp, data1.acc_y, data1.timestamp, data1.acc_z);
legend('acc.x','acc.y','acc.z');
title('Acelerómetro - L2Data1');
xlabel('Tempo [ms]');
ylabel('Aceleração [m/s^2]');

subplot(2,1,2);
plot(data1.timestamp, data1.gyro_x, data1.timestamp, data1.gyro_y, data1.timestamp, data1.gyro_z);
legend('gyro.x','gyro.y','gyro.z');
title('Giroscópio - L2Data1');
xlabel('Tempo [ms]');
ylabel('Velocidade angular [rad/s]');

saveas(gcf, 'L2Data1_Sensores.png');

acc_mean = mean([data1.acc_x, data1.acc_y, data1.acc_z]);
acc_std = std([data1.acc_x, data1.acc_y, data1.acc_z]);
acc_var = var([data1.acc_x, data1.acc_y, data1.acc_z]);

gyro_mean = mean([data1.gyro_x, data1.gyro_y, data1.gyro_z]);
gyro_std = std([data1.gyro_x, data1.gyro_y, data1.gyro_z]);
gyro_var = var([data1.gyro_x, data1.gyro_y, data1.gyro_z]);

disp('Acelerómetro (L2Data1)');
disp(acc_mean); disp(acc_std); disp(acc_var);

disp('Giroscópio (L2Data1)');
disp(gyro_mean); disp(gyro_std); disp(gyro_var);
