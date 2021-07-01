%%% Computation of Integration of Dead Reckonging and GNSS %%%
%%% Define the position at all epochs from all data %%%
function Integration_Result = Calculation_Integration
Define_Constants  % Import 'Define Constants m file and this is useful to calculate the results
Dead_reckoning = csvread('Data_File\Dead_reckoning.csv');  % Import data from the Dead reckoning csv file
[i, ~] = size(Dead_reckoning); % i means epoch

% Explain the format Integration_Result
% rows = epoch, which means time step
% 1st column = Time (s)
% 2nd column = Latitude (degrees)
% 3rd column = Longitude (degrees)
% 4th column = Velocity_x (North) (m/s)
% 5th column = Velocity_y (East) (m/s)
% 6th column = Heading (degrees)
Integration_Result = zeros(i, 6); % row = each epoch (times) and column = same format of the example output profile 

% Define the initial parameters from the Calculation_GNSS function
initial_para_GNSS = Calculation_GNSS_with_Outlier_detection;

% Define the initial parameters from the Calculation_Dead_Reckoning function
initial_para_Dead_Reckoning = Calculation_Dead_Reckoning;

% Store time in column 1
Integration_Result(:, 1) = initial_para_GNSS(:, 1);
Integration_Result(:, 6) = initial_para_Dead_Reckoning(:, 6);
% First, define several parameters
S_DR = 0.01; % velocity error variance with PSD
tau_s = 0.5; % propagation interval = 0.5
sigma_v = 0.1; % initial velocity uncertainty is 0.1m/s in each direction
sigma_r = 10; % initial position uncertainty is 10m per direction
sigma_Gr = 10; % GNSS poistion noise std on all pseudo-range measurements 
sigma_Gv = 0.05; % GNSS velocity noise std on all on all pseudo-range rate measurements 
height = initial_para_GNSS(:,7); % height computed from GNSS
Latitude = initial_para_Dead_Reckoning(:,2)*deg_to_rad; % latitude computed from Dead reckoning
% The radii of curvature may be computed from the latitude using the Matlab
% function Radii_of_curvature
[R_N_initial,R_E_initial]= Radii_of_curvature(Latitude(1));
P_k_minus1 = [sigma_v^2 0 0 0;...
             0 sigma_v^2 0 0;...
             0 0 sigma_r^2/(R_N_initial+height(1))^2 0;... 
             0 0 0 sigma_r^2/(((R_E_initial+height(1))^2)*cos(Latitude(1))^2)];
% implement a 4-state Kalman filter estimating north and east DR velocity error, 
% DR latitude error and DR longitude error
x_k_minus1 = [0;0;0;0];
% Initialize the integration result with Dead Reckoning result
Integration_Result(1,:) = initial_para_Dead_Reckoning(1,:); 
% Then follow the ten steps of the Kalman filter as follows:
for k=2:i
    % Compute the transition matrix
    [R_N,R_E]= Radii_of_curvature(Latitude(k));
    phi_k_minus1 = [1 0 0 0;...
                    0 1 0 0;...
                    tau_s/(R_N+height(k-1)) 0 1 0;...
                    0 tau_s/((R_E+height(k-1))*cos(Latitude(k-1))) 0 1];
    % Compute the system noise covariance matrix
    Q_k_minus1 = [S_DR*tau_s 0 1/2*S_DR*tau_s^2/(R_N+height(k-1)) 0;...
                  0 S_DR*tau_s 0 1/2*S_DR*tau_s^2/((R_E+height(k-1))*cos(Latitude(k-1)));...
                  1/2*S_DR*tau_s^2/(R_N+height(k-1)) 0 1/3*S_DR*tau_s^3/(R_N+height(k-1))^2 0;...
                  0 1/2*S_DR*tau_s^2/((R_E+height(k-1))*cos(Latitude(k-1))) 0 ...
                  1/3*S_DR*tau_s^3/((R_E+height(k-1))^2*cos(Latitude(k-1))^2)];
    % Propagate the state estimates
    x_k_minus = phi_k_minus1*x_k_minus1;
    % Propagate the error covariance matrix
    P_k_minus = phi_k_minus1*P_k_minus1*phi_k_minus1'+Q_k_minus1;
    % Compute the measurement matrix
    H_k = [0 0 -1 0;...
           0 0 0 -1;...
           -1 0 0 0;...
           0 -1 0 0];
    % Compute the measurement noise covariance matrix
    R_k = [sigma_Gr^2/(R_N+height(k))^2 0 0 0;...
           0 sigma_Gr^2/((R_E+height(k))^2*cos(Latitude(k))^2) 0 0;...
           0 0 sigma_Gv^2 0;...
           0 0 0 sigma_Gv^2];
    % Compute the Kalman gain matrix
    K_k = P_k_minus*H_k'/(H_k*P_k_minus*H_k'+R_k); 
    % Formulate the measurement innovation vector
    delta_z_k_minus = [initial_para_GNSS(k,2)*deg_to_rad - initial_para_Dead_Reckoning(k,2)*deg_to_rad;...
                       initial_para_GNSS(k,3)*deg_to_rad - initial_para_Dead_Reckoning(k,3)*deg_to_rad;...
                       initial_para_GNSS(k,4) - initial_para_Dead_Reckoning(k,4);...
                       initial_para_GNSS(k,5) - initial_para_Dead_Reckoning(k,5)]...
                     - H_k*x_k_minus;
    % Update the state estimates
    x_k_plus = x_k_minus+K_k*delta_z_k_minus;
    % Update the error covariance matrix
    P_k_plus = (eye(4,4) - K_k*H_k)*P_k_minus;
    
    % store the value in integration result
    % Geodetic latitude
    Integration_Result(k,2) = initial_para_Dead_Reckoning(k,2)*deg_to_rad - x_k_plus(3);
    % Geodetic longitude 
    Integration_Result(k,3) = initial_para_Dead_Reckoning(k,3)*deg_to_rad - x_k_plus(4);
    % North velocity
    Integration_Result(k,4) = initial_para_Dead_Reckoning(k,4) - x_k_plus(1); 
    % East velocity
    Integration_Result(k,5) = initial_para_Dead_Reckoning(k,5) - x_k_plus(2); 
    
    % update parameters
    P_k_minus1 = P_k_plus;
    x_k_minus1 = x_k_plus;
end
Integration_Result(2:end,2) = Integration_Result(2:end,2)*rad_to_deg;
Integration_Result(2:end,3) = Integration_Result(2:end,3)*rad_to_deg;
% csvwrite can only store value in short format, so we use dlmwrite
%dlmwrite('Integration_computation.csv', Integration_Result, 'delimiter', ',', 'precision', 9);
end