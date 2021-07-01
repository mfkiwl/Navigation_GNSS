%% Computation of GNSS and Kalman filter %%
%% Define the position at all epochs from all data %%
function GNSS_Results = Calculation_GNSS_with_Outlier_detection
    %% Define position and velocity from pseudo ranges data
    Pseudo_ranges = load('Data_File/Pseudo_ranges.csv'); % Import data from the Pesudo range csv file
    Pseudo_range_rate = load('Data_File/Pseudo_range_rates.csv'); % Import data from the Pesudo range rates csv file
    Define_Constants; % Import 'Define Constants m file and this is useful to calculate the results'

    % Define the parameters (current_satellite, num_of_satellites)
    current_satellite = Pseudo_ranges(1, 2:end); % Means each satellite number
    num_of_satellites = size(current_satellite, 2); % Means the number of satellite

    % Make an array to save the position and velocity of the satellites at 0 epoch
    pos_vel_satellite = zeros(3, num_of_satellites, 2);
    for i = 1:num_of_satellites
        inform_satellite = current_satellite(1, i); % Define the current satellite number
        % Using the Satellite_position_and_velocity.m file to define the position and velocity of the current satellite
        [Pos_satellite, Vel_satellite] = Satellite_position_and_velocity(0, inform_satellite);
        pos_vel_satellite(:, i, 1) = Pos_satellite; % Store the position of the current satellite into the array
        pos_vel_satellite(:, i, 2) = Vel_satellite; % Store the velocity of the current satellite into the array
    end

    % From now, we need to define the remaining of the columns what meaning is
    % Define the initial state vector estimate(x_0) and error covariance matrix(P)
    [state_vector_est_x, e_covariance_m_P] = Initial_GNSS_Position(Pseudo_ranges(2, 2:end), Pseudo_range_rate(2, 2:end), pos_vel_satellite);

    % Let's Define the array of the GNSS_Results format 
    % Define the number of epoch
    num_of_epoch = size(Pseudo_ranges, 1); % Same as the number of times
    GNSS_Results = zeros(num_of_epoch - 1, 7); % Design (num_of_epoch x 7) matrix

    %% Save all computation results of the satellites at the first epoch %%
    %   L_b           latitude (rad)
    %   lambda_b      longitude (rad)
    %   h_b           height (m)
    %   v_eb_n        velocity of body frame w.r.t. ECEF frame, resolved along north, east, and down (m/s) 3x1 column vector
    %   x_plus_k      state estimate vector (8x1)matrix [position-x,y,z; velocity-x,y,z; clock-offset; clock-drift] 
    [L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(state_vector_est_x(1:3), state_vector_est_x(4:6));
    % The unit of 1st and 2nd rows is radian, so it should convert from radian to degree by using the Define_Constant.m file
    Latitude = L_b*rad_to_deg; % Latitude (rad -> degrees)
    Longtitude = lambda_b*rad_to_deg; % Longitude (rad -> degrees
    Heights = h_b; % Height (m)
    GNSS_Results(1, 2:5) = [Latitude,Longtitude, v_eb_n(1:2)'];
    GNSS_Results(1, 7) = Heights;
    
    % Initialised and save the first position and velocity of the satellites in the array format
    % GNSS_Results(1, 2:5) = [Latitude, Longitude, v_eb_n(1:2)']; 
    % GNSS_Results(1, 7) = Height;
    % Explain the format of Final_Solutions          % Explain the format GNSS_Result (same as Example_Output_Profile)
    % rows = epoch, which means each time step       % rows = epoch, which means time step
    % 1st column = Time (s)                          % 1st column = Time (s)
    % 2nd column = Latitude (degrees)                % 2nd column = Latitude (degrees)
    % 3rd column = Longitude (degrees)               % 3rd column = Longitude (degrees)
    % 4th column = Heading (degrees)                 % 4th column = Velocity_x (North) (m/s)
    % 5th column = Velocity_x (North) (m/s)          % 5th column = Velocity_y (East) (m/s)
    % 6th column = Velocity_y (East) (m/s)           % 6th column = Heading (degrees)
    % 7th column = Velocity_z (m/s)                  % 7th column = Height (m) 
    % Reminder ! -> 7th column Height doesn't need to show, but we just save it
    % because we want to check the value is reasonable or not and also the height will use to compute a Dead Reckoning and Sensor integration later
    
    % From now, define the measured pesudo ranges and pesudo range rates for every epoch
    for i = 1:num_of_epoch - 1
        times = Pseudo_ranges(i + 1, 1); % Define the times exactly
        Pseudo_ranges_data = Pseudo_ranges(i + 1, 2:end); % Define measured pseudo ranges
        Pseudo_range_rates_data = Pseudo_range_rate(i + 1, 2:end); % Define measured pseudo range rates   

         % State the outputs for using the above all parameters
        % Compute the updating state estimate vector and error covariance matrix for using Kalman Filter method
        [x_plus_k, P_plus_k, Check_Outlier] = Kalman_Filter(state_vector_est_x, e_covariance_m_P,...
                                                            Pseudo_ranges_data, Pseudo_range_rates_data,...
                                                            current_satellite, times);

        %%% Outlier detection should be considered in this computation %%%
        % If the computation detects the outliers, recalculate the position at the current epoch without the measurement that had the largest residual,
        % retaining any outlying measurements with smaller residuals. After that, should repeat the outlier detection test and keep going on the computation for remaining epoch
        if ~isempty(Check_Outlier)
            % Check the outlier has contaminated the position solution or not
            [New_pseudo_ranges_data,New_pseudo_range_rates_data, New_satellite_number] = Outlier_Detections(Pseudo_ranges_data, Pseudo_range_rates_data,...
                                                                                                            current_satellite, Check_Outlier);
            % After checking the outlier detection solutions, recalculate all measurements now past the test
            [x_plus_k, P_plus_k, ~] = Kalman_Filter(state_vector_est_x, e_covariance_m_P,...
                                                    New_pseudo_ranges_data, New_pseudo_range_rates_data,...
                                                    New_satellite_number, times);
        end
        
        % Already mentioned these parameters above and the parameters should define again because the computation is repeated because of outlier detection
        [L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(x_plus_k(1:3), x_plus_k(4:6));
        % The unit of 1st and 2nd rows is radian, so it should convert from radian to degree by using the Define_Constant.m file
        Latitude = L_b * rad_to_deg; % Latitude (rad -> degrees)
        Longtitude = lambda_b * rad_to_deg; % Longitude (rad -> degrees)
        Heights = h_b; % Height (m)
        
        % Re-initialised and save the first position and velocity of the satellites in the array format
        GNSS_Results(i, 1) = times;
        GNSS_Results(i,2:5) = [Latitude, Longtitude, v_eb_n(1:2)'];

        % 1st(times), 2nd(latitude), 3rd(longitude), 4th(velocity-North), 5th(velocity-East) and 7th(Height) columns are already defined
        % Therefore, remaining column is 6th, which is a heading value. Form now, compute the 6th column(Heading)
        Heading = atan(GNSS_Results(i, 5)/ GNSS_Results(i, 4)) .* rad_to_deg;
        if GNSS_Results(i, 4) < 0  % Checking the degrees
            Heading = Heading + 180;
        end
        
        GNSS_Results(i, 6) = Heading; % Save the heading values in the 6th column
        GNSS_Results(i, 7) = Heights; % Save the height values in the 7th column
        % Final step to update the results of Kalman Filter for computing the next epoch
        state_vector_est_x = x_plus_k;
        e_covariance_m_P = P_plus_k;
    end
end

%% Initialise the Kalman Filter, this computation is followed by Workshop 1 and example of least sqaure approximation %%
%% State estimate vector and error covariance matrix %%
function [update_state_vector_est_x, update_e_covariance_m_P] = Initial_GNSS_Position(Pseudo_r, Pseudo_r_r, Pos_Vel_s)
    % Compute the linear least square approximation method to initialised the user position
    % by using the initial position of the satellite
    Define_Constants; % Import 'Define Constants m file and this is useful to calculate the results'
    position_satellite = Pos_Vel_s(:, :, 1); % Define the position of satellite by using Satellite_position_and_velocity.m file
    num_of_satellite  = size(position_satellite, 2); % Define the number of satellite
    state_vector_est_x = zeros(8, 1); % Define 8x1 matrix for state estimate vector
    new_position = zeros(size(position_satellite)); % Define the new position for using Sagnac effect compensation matrix
    
    % Compute the new position to check the correct position by using the Sagnac effect compensation matrix
    for k = 1:num_of_satellite
        r_aj = Pseudo_r(k); % Range from the approximation the user position to each satellite
        current_position = position_satellite(:, k); % Define the position of the satellites
        % Compute the Sagnac effect compensation matrix
        C_I_e = [               1    omega_ie*r_aj/c     0;...
                 -omega_ie*r_aj/c                  1     0;...
                                0                  0     1];
        new_position(:, k) = C_I_e * current_position;
    end
    % After finishing the computation of the Sagnac effect, then compute the new range
    new_ranges = Pseudo_r(1);
    Total_ranges = Pseudo_r;
    Pseudo_r(1) = [];
    new_ranges = repmat(new_ranges, [1, num_of_satellite - 1]);
    
    % Also compute the new position
    new_positions = new_position(:, 1);
    new_positions = repmat(new_positions, [1, num_of_satellite - 1]);
    new_position(:, 1) = [];
    
    % Should check the distance between the current satellite of the position and next satellite of the position 
    element_1 = (new_position - new_positions)';
    cs_ns_distance = sum((new_position(:, :) - new_positions(:, :)).^2, 1);
    element_2 = 1/2.*(new_ranges.^2 - Pseudo_r.^2 + cs_ns_distance).';
    A = element_1 \ element_2;
    state_vector_est_x(1:3) = A + new_positions(:, 1);

    % From now, we finished to initialised the position of the state vector [Position-x, Position-y, Position-z]
    % Then, define the velocity of the state vector [Velocity-North, Velocity-East, Velocity-Down]
    state_vector_est_x(4:6) = zeros(3, 1);
    % Define the clock offset and clock drift of the state estimate vector
    % Uses least-squares to estimate the receiver clock errors before initialise an extended Kalman filter
    clock_offset = 100000; % The clock offset (m) would be computed to extend the Kalman filter
    clock_drift = 200; % The clock drift (m/s) would be computed to extend the Kalman filter
    state_vector_est_x(7:8) = [clock_offset; clock_drift];
    
    % Initialise the state estimate vector when the above computation is finished at the first epoch cycle
    update_state_vector_est_x = Calculation_GNSS_Receiver_Clock(state_vector_est_x, Total_ranges, Pos_Vel_s, Pseudo_r_r);

    % Initialise the error covariance matrix at the first epoch cycle and From the coursework 1 instruction,
    sigma_pseudo_r = 10; % The noise standard deviation assumes 10 m on all pseudo range measurements
    sigma_pseudo_r_r = 0.05; % The noise standard deviation assumes 0.05 m/s on all pseudo range rates measurements
    update_e_covariance_m_P = [sigma_pseudo_r^2 * eye(3)                                        zeros(3, 5);...
                                                zeros(3)   sigma_pseudo_r_r^2 * eye(3)          zeros(3, 2);...
                                             zeros(1, 6)              sigma_pseudo_r^2                    0;...
                                             zeros(1, 7)                                 sigma_pseudo_r_r^2];    
%     update_e_covariance_m_P = [100     0     0       0      0      0      0      0;    for position-x
%                                  0   100     0       0      0      0      0      0;    for position-y
%                                  0     0   100       0      0      0      0      0;    for position-z 
%                                  0     0     0  0.0025      0      0      0      0;    for velocity-x 
%                                  0     0     0       0 0.0025      0      0      0;    for velocity-y
%                                  0     0     0       0      0 0.0025      0      0;    for velocity-z
%                                  0     0     0       0      0      0    100      0;    for clock offset
%                                  0     0     0       0      0      0      0 0.0025]    for clock drift    
end

%% Initialise the Kalman Filter, this computation is followed by Workshop 1 %%
%% State clock offset, clock drift and compute the update velocity %%     
function update_state_vector_est_x = Calculation_GNSS_Receiver_Clock(state_vector_est_x, Total_ranges, Pos_Vel_s, Pseudo_r_r)
    % Define the parameters to compute the results
    position_satellite = Pos_Vel_s(:, :, 1); % Define the position of satellite by using Satellite_position_and_velocity.m file
    velocity_satellite = Pos_Vel_s(:, :, 2); % Define the velocity of satellite by using Satellite_position_and_velocity.m file
    num_of_satellite = size(Total_ranges, 2); % Define the number of satellite
    Define_Constants; % Import 'Define Constants m file and this is useful to calculate the results

    % Check the results to converge to constant values during the computation
    while true
        %%%%%%%%%%% With respect to position and clock offset %%%%%%%%%
        % Compute the clock offset by using the predicted ranges for each satellite
        % Make an array of predicted range format
        predict_ranges = zeros(1, num_of_satellite);
        % Make a for loop and compute the predicted ranges for each statellite 
        for k = 1:num_of_satellite
            current_position = position_satellite(:, k);
            distance = current_position - state_vector_est_x(1:3);
            % Predict the ranges from the approximatie user position to each satellite
            r_aj = sqrt(distance.' * distance);
            % Compute the Sagnac effect compensation matrix
            C_I_e = [                   1    omega_ie * r_aj / c    0;...
                      -omega_ie * r_aj /c                      1    0;...
                                        0                      0    1];
            % Re-compute the predicted range
            new_distance = C_I_e * current_position - state_vector_est_x(1:3);
            predict_ranges(:, k) = sqrt(new_distance.' * new_distance);
        end

        % After finishing the computation of the predict ranges from the approximate user position to each satellite,
        % then need to compute the line-of-sight unit vector from the approximate user position to each satellite
        u_e_aj = (position_satellite - repmat(state_vector_est_x(1:3), [1, num_of_satellite]))./repmat(predict_ranges, [3,1]);
        % easy method for using approximation, but this result slightly different results 
        % so we are not using this equation. However we write this equation to help understanding 
        % u_e_aj = (current_position - state_vector_est_x(1:3)) / position_satellite;
        
        % Formulate the predicted state vector, measurement innovation vector (delta_z), and measurement matrix (H_e_G)
        % Compute the measurement innovation vector (delta_z)
        % Input = Total_ranges: the measured pseudo range from each satellite to the user antenna  
        %         predict_ranges: Predict the ranges from the approximate user position to each satellite
        %         delta_rho_a_c: Predicted receiver clock offset
        delta_rho_a_c = ones(1, num_of_satellite) .* state_vector_est_x(4);
        delta_z = Total_ranges - predict_ranges - delta_rho_a_c;
        
        % Compute the measurement matrix (H_e_G)
        H_e_G = [-u_e_aj.'    ones(num_of_satellite, 1)]; % Concatenates the line-of-sight matrix with ones(8x1) matrix
        
        % After finishing the compuation of the measurement innovation vector (delta_z) and matrix (H_e_G),
        % then compute the position and receiver clock offset using unweighted least sqaures
        x_plus_pos = state_vector_est_x(1:4) + (H_e_G.' * H_e_G)^-1 * H_e_G.' * delta_z.';

        %%%%%%%%%%% With respect to velocity and clock drift %%%%%%%%%
        % Compute the clock drift by using the predicted ranges for each satellite
        % Make an array of predicted range rates format
        predict_range_rates= zeros(1, num_of_satellite);
        % Make a for loop and compute the predicted ranges for each statellite 
        for k = 1:num_of_satellite
            r_aj = Total_ranges(1, k);
            current_position = position_satellite(:, k);
            current_velocity = velocity_satellite(:, k);
            
            % Compute the predicted range rates (r_aj_dot)
            % Following the equation (9) in Workshop 1,
            element_1 = current_velocity + Omega_ie * current_position;
            element_2 = state_vector_est_x(5:7) + Omega_ie * state_vector_est_x(1:3);
            u_e_aj_vel = u_e_aj(:, k).';
            % Compute the Sagnac effect compensation matrix for the velocity
            C_I_e = [                   1    omega_ie * r_aj / c   0;...
                     -omega_ie * r_aj / c                      1   0;...
                                        0                      0   1];
            predict_range_rates(1, k) = u_e_aj_vel * (C_I_e * element_1 - element_2);
        end

        % Formulate the predicted state vector, measurement innovation vector (delta_z), and measurement matrix (H_e_G)
        % Compute the measurement innovation vector (delta_z)
        % Input = Total_ranges: the measured pseudo range from each satellite to the user antenna  
        %         predict_ranges: Predict the ranges from the approximate user position to each satellite
        %         delta_rho_a_c: Predicted receiver clock offset
        rho_j_a_dot = Pseudo_r_r;
        delta_rho_a_c_dot = ones(1, num_of_satellite) * state_vector_est_x(8);
        delta_z_vel = rho_j_a_dot - predict_range_rates - delta_rho_a_c_dot;
        
        % Compute the measurement matrix (H_e_G) for velocity
        H_e_G_vel = [-u_e_aj.'    ones(num_of_satellite, 1)] ; % Concatenates the line-of-sight matrix with ones(8x1) matrix
        
        % Update the state estimate vector
        x_plus_vel = state_vector_est_x(5:8) + (H_e_G_vel.' * H_e_G_vel)^-1 * H_e_G_vel.' * delta_z_vel.';

        % Update the state estimate vector
        new_state_vector_est_x = [x_plus_pos;...
                                  x_plus_vel];  % Make a 8x1 vector

        % Check the difference bewteen old version of the state vector and new version of the state vector
        %converge_test = sqrt(sum((state_vector_est_x - new_state_vector_est_x).^2, 1));
        converge_test = norm(state_vector_est_x - new_state_vector_est_x);
        % If the value is bigger than 0.0001, keep runing the while loop to make convergence until less than 0.0001
        % If the value is less than 0,0001, stop and update the new state estimate vector  
        if converge_test >= 0.01  % Reduce the error
            state_vector_est_x = new_state_vector_est_x;
        else
            break;
        end
    end
    % Final update and use these vector (8x1)
    update_state_vector_est_x = [x_plus_pos;...
                                 x_plus_vel];
end

%% We finished to compute the Initialised Kalman Filter at the 1st epoch %%
%% Now the Kalman Filter at all epochs should be computed to get position and velocity %%
function [x_plus_k, P_plus_k, Outlier_Detection] = Kalman_Filter(state_vector_est_x, e_covariance_m_P,...
                                                                 New_pseudo_ranges_data, New_pseudo_range_rates_data,...
                                                                 New_satellite_number, times)
    % Define the parameters to compute the results
    Define_Constants; % Import 'Define Constants m file and this is useful to calculate the results'
    
    % Compute the every epoch (time) and Following the workshop 2 Task 2A
    %% Kalman Filter Step 1 %%
    % Compute the transition matrix (input: Identity matrix and tau_s)
    tau_s = 0.5; % Propagation interval tau_s(s)
    transition_m = [    eye(3) tau_s * eye(3)  zeros(3, 1)  zeros(3, 1);...
                      zeros(3)         eye(3)  zeros(3, 1)  zeros(3, 1);...
                   zeros(1, 3)    zeros(1, 3)            1       tau_s;...
                   zeros(1, 3)    zeros(1, 3)            0           1];

    %% Kalman Filter Step 2 %%
    % Compute the system noise covariance matrix 
    % (input: the acceleration power spectral density (PSD) S_e_a = 3 (m^2/s^3)
    S_e_a = 5 ; % Value of the acceleration for a pedestrian as workshop 2 and this accounts for the propagation of the system noise onto the position state during the propagation interval
    S_a_cphi = 0.01; % Clock phase PSD for the GNSS receiver clock
    S_a_cf =  0.04; % Clock frequency PSD for the GNSS receiver clock
    Q_k = [(S_e_a * tau_s^3 * eye(3))/3  (S_e_a * tau_s^2 * eye(3))/2                               zeros(3, 1)           zeros(3, 1);...
           (S_e_a * tau_s^2 * eye(3))/2      (S_e_a * tau_s * eye(3))                               zeros(3, 1)           zeros(3, 1);...
                            zeros(1, 3)                   zeros(1, 3) (S_a_cphi * tau_s + (S_a_cf * tau_s^3)/3)  (S_a_cf * tau_s^2)/2;...
                            zeros(1, 3)                   zeros(1, 3)                      (S_a_cf * tau_s^2)/2        S_a_cf * tau_s];
                        
    %% Kalman Filter Step 3 %%
    % Propagate the state estimate
    x_minus_k = transition_m * state_vector_est_x;
    
    %% Kalman Filter Step 4 %%
    % Propagate the error covariance matrix
    P_minus_k =  transition_m * e_covariance_m_P * transition_m.' + Q_k;

    % After finishing the computation of the state estimate vector and the error covariance matrix,
    % then need to compute the line-of-sight unit vector from the approximate user position to each satellite
    % (Following the Workshop 2 Task 2A - f) and g))
    num_of_satellite = size(New_satellite_number, 2);
    Pos_Vel_s = zeros(3, num_of_satellite, 2);
    % Compute the position and velocity of each satellite at the current epoch
    for k = 1:num_of_satellite
        current_satellite = New_satellite_number(1, k); % Define the current satellite
        %   position_satellite       ECEF satellite position (m) vector
        %   sat_v_es_e               ECEF satellite velocity (m/s) vector
        [position_satellite, sat_v_es_e] = Satellite_position_and_velocity(times, current_satellite);
        Pos_Vel_s(:, k, 1) = position_satellite;
        Pos_Vel_s(:, k, 2) = sat_v_es_e;
    end

    % Compute the predicted range at the current epoch
    predict_ranges = zeros(1, num_of_satellite);
    for k = 1:num_of_satellite
        % Initialised the predicted range
        current_position = Pos_Vel_s(:, k, 1);
        distance = current_position - x_minus_k(1:3, :);
        r_aj = sqrt(distance.' * distance);

        % Compute the Sagnac effect compensation matrix
        C_I_e = [                  1      omega_ie * r_aj/c    0;...
                  -omega_ie * r_aj/c                      1    0;...
                                   0                      0    1];
        % Re-compute the predicted range
        new_distance = C_I_e * current_position - x_minus_k(1:3, :);
        predict_ranges(:, k) = sqrt(new_distance.' * new_distance);
    end
    % From now, compute the line of sight unit vector from the approximate user position to each satellite
    u_e_aj = (Pos_Vel_s(:, :, 1) - repmat(x_minus_k(1:3, :), [1, num_of_satellite]))./repmat(predict_ranges, [3, 1]);

    %% Kalman Filter Step 5 %%
    % Compute the measurement matrix
    H_k = [                 -u_e_aj.'  zeros(num_of_satellite, 3)   ones(num_of_satellite, 1)  zeros(num_of_satellite, 1);...
           zeros(num_of_satellite, 3)                   -u_e_aj.'  zeros(num_of_satellite, 1)   ones(num_of_satellite, 1)];
       
    %% Kalman Filter Step 6 %%
    % Compute the measurement noise covariance matrix
    sigma_p = 2; % Code tracking and multipath error standard deviation for pseudo range measurement (m)
    sigma_r = 0.02; % Range rate tracking and multipath error standard deviation for pseudo range rate measurement (m/s)
    R_k = [sigma_p^2 * eye((num_of_satellite),(num_of_satellite))                zeros(num_of_satellite, num_of_satellite);...
                        zeros(num_of_satellite, num_of_satellite)     sigma_r^2 * eye(num_of_satellite, num_of_satellite)];
    R_k = R_k + [zeros(num_of_satellite, 2 * (num_of_satellite));...
                 zeros(num_of_satellite,  2 * (num_of_satellite) - 1)   ones(num_of_satellite, 1)];
    R_k(2 * (num_of_satellite), 2 * (num_of_satellite)) = 0;%R_k(2 * (num_of_satellite), 2 * (num_of_satellite)) - 1;
    % This method for computing the measurement error covariance follows
    % the COMP0130 Least Sqaure Step by Step.pdf 
    
    %% Kalman Filter Step 7 %%
    % Compute the Kalman gain matrix
    K_k = P_minus_k * H_k.' * (H_k * P_minus_k * H_k.' + R_k)^-1;

    % Compute the predict range rate for the velocity
    r_aj_dot = zeros(1, num_of_satellite);
    for k = 1:num_of_satellite
        % Define the position and velocity for each satellite for using the Satellite_position_and_velocity.m file
        position_satellite = Pos_Vel_s(:, k, 1);
        velocity_satellite = Pos_Vel_s(:, k, 2);

        % Define the propagated state vector for the position and velocity for each satellite
        r_e_ea = x_minus_k(1:3, :); % Predicted Cartesian ECEF user position
        v_e_ea = x_minus_k(4:6, :); % Predicted Cartesian ECEF user velocity

        % Predict the range rates from the approximate user position to each satellite
        r_aj = predict_ranges(1, k);

        % Define the line-of-sight unit vector from the approximate user position to each satellite
        element_1 = velocity_satellite + Omega_ie * position_satellite;
        element_2 = v_e_ea + Omega_ie * r_e_ea;
        LoS = u_e_aj(:, k).'; % LoS means the line-of-sight
        C_I_e = [                   1  omega_ie * r_aj / c    0;...
                 -omega_ie * r_aj / c                    1    0;...
                                    0                    0    1];
        r_aj_dot(1, k) = LoS * (C_I_e * element_1 - element_2);
    end

    %% Kalman Filter Step 8 %%
    % Formulate the measurement innovation vector
    delta_z = zeros(1, 2 * (num_of_satellite)); % Make an array size of the measurement innovation vector
    delta_rho_a_c_dot = ones(1, num_of_satellite) * x_minus_k(end - 1);
    delta_rho_a_c_dot_again = ones(1, num_of_satellite) * x_minus_k(end);
    delta_z(1, 1:num_of_satellite) = New_pseudo_ranges_data - predict_ranges - delta_rho_a_c_dot;
    delta_z(1, num_of_satellite + 1:end) = New_pseudo_range_rates_data - r_aj_dot - delta_rho_a_c_dot_again;

    %% Kalman Filter Step 9 %%
    % Update the state estimates
    x_plus_k = x_minus_k + K_k * delta_z.';

    %% Kalman Filter Step 10 %%
    % Update the error covariance matrix
    P_plus_k = (eye(size(K_k * H_k))- K_k * H_k) * P_minus_k;

    %% Finished the computation of Kalman Filter at the epoch %%
    %% Check the residual based outlier detection for the position at the epoch %%
    % Define inputs: measurement matrix (H_e_G_outlier), measurement innovation matrix (delta_z_outlier),...
    %                ionosphere error (sigma_iono), troposhpere error (sigma_tropo),
    %                outlier detection threshold (T)
    sigma_iono = 2; % Define the residual ionosphere error standard deviation at zenith
    sigma_tropo = 0.2; % Define the residual troposhpere error standard deviation at zenith
    T = 6; % Define the outlier detection threshold, which is based on the workshop 1 Task 3 - (c)
    H_e_G_outlier = [H_k(1:num_of_satellite , 1:3), ones(num_of_satellite, 1)] ;
    delta_z_outlier = delta_z(1, 1:num_of_satellite) ;
    
    %% Outlier Detection Step 1 %%
    % Compute the residuals vector
    v = (H_e_G_outlier * ((H_e_G_outlier.' * H_e_G_outlier)^-1) * H_e_G_outlier.' - eye(num_of_satellite)) * delta_z_outlier.';

    % Consider for the position (x,y,z) in line-of-sight unit vector, which means latitude and longtitude, respectively, so we need to measure both error at Zenith
    % The angle between the user and satellite can compute the sine function sin(z_e_as, r_as)
    SD_zenith = repmat((sigma_iono + sigma_tropo), 1, num_of_satellite); % Total residual error for the position is sum with ionosphere and troposphere
    Angle = u_e_aj(3, :)./r_aj;
    SD_measurement_error = SD_zenith./sin(Angle); % Reference by Lecture 1D Slide pg 21
    % Compute the measurement error variance = (Measurement error of standard deviation)^2
    e_variance = (SD_measurement_error).^2;

    %% Outlier Detection Step 2 %%
    % Compute the residuals covariance matrix
    C_v = eye(num_of_satellite) - H_e_G_outlier * (H_e_G_outlier.' * H_e_G_outlier)^-1 * H_e_G_outlier.';

   %% Outlier Detection Step 3 %%
    % Compute the normalised residuals and compare each with a threshold
    C_v_jj = diag(C_v)' .* e_variance; % Define the diagonal element of the residuals covariance matrix
    w_i = v' ./ sqrt(C_v_jj); % Normalised the residual, based on Least Sqaure Step-by-Step.pdf (Section: Measurement error hypothesis testing)

    %% Outlier Detection Step 4 %%
    % Check the condition when the measurement j is an outlier with Threshold
    check_condition = abs(w_i) > T;  % same as abs(v) > sqrt(C_v_jj).* T;
    outlier = abs(w_i) - T;
    if outlier > 0 % If the condition of the outlier detection is satisfied,
        Outlier_Detection = check_condition; % Save the values and keep going to compute for the next epoch
    else                        % If not, the outlier has contaminated the position solution and have the largest residual, so
        Outlier_Detection = []; % Remove the measurement and repeat the emasurement to find that all measurements
    end
end

%% Handle the outlier detection %%
%% After Computing the Kalman Filter, if the measurement had the largest residual, then make empty array to recalculate the position at the epoch %%
function [update_r, update_r_r, update_satellite] = Outlier_Detections(previous_r, previous_r_r, previous_satellite, Check_Outlier)
    % Errors in individual measurements or their measurement models can be identified from the residuals
    % Define the predicted ranges 
    update_satellite = previous_satellite;
    update_r = previous_r;
    update_r_r = previous_r_r;
    num_of_satellites = size(previous_satellite, 2);

    occur_outlier = 0; % Define the initial Outlier detection
    % Removing the measurement from the outlier detection of the satellite
    % and changes the range data to compute for solving right results
    for k = 1:num_of_satellites
        if Check_Outlier(k) == 1 % If the outlier detection occurs, then use this loop
            occur_outlier = occur_outlier + 1; % Comfirm the detection outliers on some of the other measurements
            update_satellite(k - occur_outlier) = []; % Found the outlier detection of the satellite number, which has the largest residual
            update_r_r(k - occur_outlier) = []; % Save new range rate values and remove the measurement and repeat the emasurement to find a new range rate
            update_r(k - occur_outlier) = []; % Remove the measurement and repeat the mesurement to find a new range
        end
    end

end
