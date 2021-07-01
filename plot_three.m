initial_para_GNSS = Calculation_GNSS_with_Outlier_detection;
initial_para_Dead_Reckoning = Calculation_Dead_Reckoning;
initial_para_Integration = Calculation_Integration;

% time parameter
time = initial_para_GNSS(:,1);
% north velocity parameters
N_velocity_GNSS = initial_para_GNSS(:, 4);
N_velocity_Dead_Reckoning = initial_para_Dead_Reckoning(:, 4);
N_velocity_Integration = initial_para_Integration(:, 4);
% east velocity parameters
E_velocity_GNSS = initial_para_GNSS(:, 5);
E_velocity_Dead_Reckoning = initial_para_Dead_Reckoning(:, 5);
E_velocity_Integration = initial_para_Integration(:, 5);
% latitude parameters
lat_GNSS = initial_para_GNSS(:, 2);
lat_Dead_Reckoning = initial_para_Dead_Reckoning(:, 2);
lat_Integration = initial_para_Integration(:, 2);
% longitude parameters
long_GNSS = initial_para_GNSS(:, 3);
long_Dead_Reckoning = initial_para_Dead_Reckoning(:, 3);
long_Integration = initial_para_Integration(:, 3);
% north velocity
figure,
plot(time,N_velocity_GNSS,'r','LineWidth',5);

hold on;

plot(time,N_velocity_Dead_Reckoning,'m');

hold on;

plot(time,N_velocity_Integration,'b');

xlabel('Time (seconds)')
ylabel('Velocity North (m/s)')
legend({'GNSS','Dead Reckoning','Integration'}, 'FontSize',8)
title('North Velocity Solutions')
% east velocity
figure,
plot(time,E_velocity_GNSS,'r','LineWidth',1.5);

hold on;

plot(time,E_velocity_Dead_Reckoning,'m');

hold on;

plot(time,E_velocity_Integration,'b');

xlabel('Time (seconds)')
ylabel('Velocity Easy (m/s)')
legend({'GNSS','Dead Reckoning','Integration'}, 'FontSize',8)
title('East Velocity Solutions')
% latitude
figure,
plot(time,lat_GNSS,'r');

hold on;

plot(time,lat_Dead_Reckoning,'m');

hold on;

plot(time,lat_Integration,'b');

xlabel('Time (seconds)')
ylabel('Degrees')
legend({'GNSS','Dead Reckoning','Integration'}, 'FontSize',8,'Location','northwest')
title('Latitude Solutions')

% longitude
figure,
plot(time,long_GNSS,'r');

hold on;

plot(time,long_Dead_Reckoning,'m');

hold on;

plot(time,long_Integration,'b');

xlabel('Time (seconds)')
ylabel('Degrees')
legend({'GNSS','Dead Reckoning','Integration'}, 'FontSize',8,'Location','southeast')
title('Longitude Solutions')
% position
figure,
plot(long_GNSS,lat_GNSS,'r');

hold on;

plot(long_Dead_Reckoning,lat_Dead_Reckoning,'m');

hold on;

plot(long_Integration,lat_Integration,'b');

xlabel('Longitude{Degrees)')
ylabel('Latitude (Degrees)')
legend({'GNSS','Dead Reckoning','Integration'}, 'FontSize',8,'Location','southeast')
title('Position Solutions')