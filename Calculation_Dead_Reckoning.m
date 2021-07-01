%%% Computation of Dead Reckonging and Integration of Gyro-Magnetometer %%%
%%% Define the position at all epochs from all data %%%
function Dead_Reckoning_Result = Calculation_Dead_Reckoning
Define_Constants  % Import 'Define Constants m file and this is useful to calculate the results
Dead_reckoning = csvread('Data_File\Dead_reckoning.csv');  % Import data from the Dead reckoning csv file
[i, ~] = size(Dead_reckoning); % i means epoch
% Make array format to save the results after finishing the computation

% Explain the format Dead_Reckoning_Result
% rows = epoch, which means time step
% 1st column = Time (s)
% 2nd column = Latitude (degrees)
% 3rd column = Longitude (degrees)
% 4th column = Velocity_x (North) (m/s)
% 5th column = Velocity_y (East) (m/s)
% 6th column = Heading (degrees)
Dead_Reckoning_Result = zeros(i, 6); % row = each epoch (times) and column = same format of the example output profile 

% Let's define the results profile for time(s) because it's same as the data
Dead_Reckoning_Result(:, 1) = Dead_reckoning(:, 1); % 1st column = Time(s) from Dead reckoning data
time = Dead_reckoning(:, 1);
% Define the initial parameters from the Calculation_GNSS function
initial_para = Calculation_GNSS_with_Outlier_detection;
%% Define several useful parameters

% Define driving speed for the robotic lawnmower
% Since it has four wheeels and the rear wheels are the driving wheels,
% we denote the average speed of rear wheels as the driving speed
driving_speed = 0.5*(Dead_reckoning(:,4)+Dead_reckoning(:, 5));

% Define gyroscope angular rate
gyroscope_angular_rate = Dead_reckoning(:,6);

% Define heading measurements in degrees from the magnetic compass
magnetic_heading = Dead_reckoning(:,7);

%% Integrate Gyro-Magnetometer to get the corrected heading
% First, we calculate gyro-heading from angular rate

% Initalize gyro-heading
gyro_heading = zeros(i,1);
% Define initial gyro-heading: when we look at the Dead reckoning csv file we
% find that the initial gyroscope angular rate is 0. Thus, the initial
% gyro-heading is 0
gyro_heading(1) = 0;
% Find total gyro-heading on each time epoch: compute the integral of 
% angular rate in the horizontal plane.
for j = 2:i
    % time interval is 0.5 based on the Dead reckoning csv file
    gyro_heading(j) = gyro_heading(j-1)+0.5*gyroscope_angular_rate(j);
end

% Then we find the integration of Gyro-Magnetometer
% First, define several parameters
% Convert deg to rad or rad to deg using Define_Constants
bias_std = 1*deg_to_rad; % A bias standard deviation of 1 degree per second
magnet_error_std = 4*deg_to_rad; % a noise-like error with a standard deviation of 4 degrees
magnetic_heading = magnetic_heading*deg_to_rad; % convert to radian
noise_std = 10^-4;% A noise standard deviation of 10-4 radians per second
S_rg = 3*10^-6; % Gyro random noise with power spectral density
S_bgd = 3*10^-6; % Gyro bias variation with PSD
tau_s = 0.5; % propagation interval = 0.5
P_k_minus1 = [noise_std^2 0; 0 bias_std^2];
x_k_minus1 = [0;0];
% initialize our result corrected heading
corrected_heading = zeros(i,1);
% Since there is no initial gyro-heading, we store the first magnetic
% heading as our first corrected heading.
corrected_heading(1) = magnetic_heading(1); 
% Then follow the ten steps of the Kalman filter as follows:
for k=2:i
    % Compute the transition matrix
    phi_k_minus1 = [1 tau_s; 0 1];
    % Compute the system noise covariance matrix
    Q_k_minus1 = [S_rg*tau_s+1/3*S_bgd*tau_s^3 1/2*S_bgd*tau_s^2;...
        1/2*S_bgd*tau_s^2 S_bgd*tau_s];
    % Propagate the state estimates
    x_k_minus = phi_k_minus1*x_k_minus1;
    % Propagate the error covariance matrix
    P_k_minus = phi_k_minus1*P_k_minus1*phi_k_minus1'+Q_k_minus1;
    % Compute the measurement matrix
    H_k = [-1 0];
    % Compute the measurement noise covariance matrix
    R_k = [magnet_error_std^2];
    % Compute the Kalman gain matrix
    K_k = P_k_minus*H_k'/(H_k*P_k_minus*H_k'+R_k); 
    % Formulate the measurement innovation vector
    delta_z_k_minus = [magnetic_heading(k) - gyro_heading(k)] - H_k*x_k_minus;
    % Update the state estimates
    x_k_plus = x_k_minus+K_k*delta_z_k_minus;
    % Update the error covariance matrix
    P_k_plus = (eye(2) - K_k*H_k)*P_k_minus;
    
    % store the value in corrected_heading
    corrected_heading(k) = gyro_heading(k)-x_k_plus(1);
    
    % update parameters
    P_k_minus1 = P_k_plus;
    x_k_minus1 = x_k_plus;
end

% store corrected heading in degrees to Dead_Reckoning_Result
Dead_Reckoning_Result(:,6) = corrected_heading*rad_to_deg;
%% Use corrected heading to find the dead Reckoning solution
% Define geodetic height for dead reckoning
height = initial_para(:,7);
% Define initial position from the Calculation_GNSS function
Dead_Reckoning_Result(1,2) = initial_para(1,2)*deg_to_rad; % initial Geodetic latitude 
Dead_Reckoning_Result(1,3) = initial_para(1,3)*deg_to_rad; % initial Geodetic longitude
% Define initial velocity: when we look at the Dead reckoning csv file we 
% find that the initial speed is 0. Thus, the initial velocity is 0.
Dead_Reckoning_Result(1,4) = 0; % initial North velocity
Dead_Reckoning_Result(1,5) = 0; % initial East velocity
% Find geodetic latitude, longitude and North and East velocity
for k = 2:i
    % The average velocity between epochs k−1 and k
    v_N_k = 1/2*(cos(corrected_heading(k))+cos(corrected_heading(k-1)))*...
        driving_speed(k);
    v_E_k = 1/2*(sin(corrected_heading(k))+sin(corrected_heading(k-1)))*...
        driving_speed(k);
    
    % The radii of curvature may be computed from the latitude using the Matlab
    % function Radii_of_curvature
    [R_N,R_E]= Radii_of_curvature(Dead_Reckoning_Result(k-1,2));
    % Update the latitude, L_k(Dead_Reckoning_Result(:,2)), 
    % and tue longitude, lambda_k(Dead_Reckoning_Result(:,3))
    Dead_Reckoning_Result(k,2) = Dead_Reckoning_Result(k-1,2) + (v_N_k*0.5)...
       /(R_N+height(k));
    Dead_Reckoning_Result(k,3) = Dead_Reckoning_Result(k-1,3) + (v_E_k*0.5)...
       /((R_E+height(k))*cos(Dead_Reckoning_Result(k,2)));
    % Update the North velocity(Dead_Reckoning_Result(:,4)) and 
    % East velocity(Dead_Reckoning_Result(:,5)) without damping
    Dead_Reckoning_Result(k,4) = 2*v_N_k - Dead_Reckoning_Result(k-1,4);
    Dead_Reckoning_Result(k,5) = 2*v_E_k - Dead_Reckoning_Result(k-1,5);
end
% convert latitude and longitude to degree
Dead_Reckoning_Result(:,2) = Dead_Reckoning_Result(:,2)*rad_to_deg;
Dead_Reckoning_Result(:,3) = Dead_Reckoning_Result(:,3)*rad_to_deg;
%csvwrite('Dead_Reckoning_computation.csv',Dead_Reckoning_Result)

end