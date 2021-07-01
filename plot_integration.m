initial_para_Integration = Calculation_Integration;
% Define parameters
time = initial_para_Integration(:,1);
N_velocity_Integration = initial_para_Integration(:, 4);
E_velocity_Integration = initial_para_Integration(:, 5);
lat_Integration = initial_para_Integration(:, 2);
long_Integration = initial_para_Integration(:, 3);
heading_Integration = initial_para_Integration(:, 6);
% plot north velocity 
figure,
plot(time, N_velocity_Integration, 'r');
xlabel('Time (seconds)')
ylabel('Velocity North (m/s)')
title('North Velocity Solutions')
% plot east velocity 
figure,
plot(time, E_velocity_Integration, 'r');
xlabel('Time (seconds)')
ylabel('Velocity North (m/s)')
title('East Velocity Solutions')
% plot geodetic latitude
figure,
plot(time, lat_Integration, 'r');
xlabel('Time (seconds)')
ylabel('Latitude (Degrees)')
title('Geodetic Latitude Solutions')
% plot geodetic longitude
figure,
plot(time,long_Integration, 'r');
xlabel('Time (seconds)')
ylabel('Longitude (Degrees)')
title('Geodetic Longitude Solutions')
% plot position
figure,
plot(long_Integration, lat_Integration, 'r');
xlabel('Longitude (Degree)')
ylabel('Latitude (Degree)')
title('Position Solutions')
legend('Trajectory of the lawnmower');
% plot heading
figure,
plot(time, heading_Integration, 'r');
xlabel('Time (seconds)')
ylabel('Heading (Degree)')
title('Heading Solutions')