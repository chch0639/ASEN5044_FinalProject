%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Mitchell Smith
% Written:  14 Dec 2017
% Revised:  14 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  ASEN 5044 - Statistical Estimation for Dynamical Systems Final
%           Project. The Kalman Filter equations were obtained from Simon, 
%           page ## (eq ##)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   xhat:   State Estimate
% 
%           k:      Discreet time Step k that EKF is on
%
%           dt:     time step EKF uses
%
%           Q:      Process noise Covariance
%
%           R:      Measurement noise Covariance
%
% Outputs:  yhat:   Random measurement given state and input
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [yhat, Hnom] = measure(xhat, u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [yhat] = noisy_measure(xhat, k, dt, Q, R)

stations = 6;           % Number of tracking stations
RE = 6378;              % Radius of Earth [km]
omegaE = 2*pi/86400;    % Rotation rate of Earth [rad/s]

% Extract States to make equations shorter and easier to read
X = xhat(1);
Xdot = xhat(2);
Y = xhat(3);
Ydot = xhat(4);

% Initialize variables in loop
yhat = zeros(3,1); 
Hnom = zeros(3,4);
Xs = zeros(stations,1); 
Xs_dot = zeros(stations,1);
Ys = zeros(stations,1); 
Ys_dot = zeros(stations,1); 
theta = zeros(stations,1);
phi = zeros(stations,1);

% Get state of each of the stations
for ii = 1:stations
    theta0 = (ii-1)*pi/3; % Initial position of each station

    % X and Y position and velocity of each staion
    Xs(ii) = RE*cos(omegaE*(k-1)*dt + theta0);
    Xs_dot(ii) = -RE*omegaE*sin(omegaE*(k-1)*dt + theta0);
    Ys(ii) = RE*sin(omegaE*(k-1)*dt + theta0);
    Ys_dot(ii) = RE*omegaE*cos(omegaE*(k-1)*dt + theta0);
    xs = [Xs(ii); Xs_dot(ii); Ys(ii); Ys_dot(ii)];
    
    
    % nonlinear measurement
    theta(ii) = atan2(Ys(ii),Xs(ii));
    phi(ii) = atan2((Y-Ys(ii)), (X-Xs(ii)));
    
    check1 = theta(ii) - pi/2 <= phi(ii) && phi(ii) <= theta(ii) + pi/2;
    check2 = theta(ii) - 5*pi/2 <= phi(ii) && phi(ii) <= theta(ii) - 3*pi/2;
    check3 = theta(ii) + 3*pi/2 <= phi(ii) && phi(ii) <= theta(ii) + 5*pi/2;
    
    if check1 || check2 || check3
        % Take measurement if satellite is in view of station
        rho = sqrt((X-Xs(ii))^2 + (Y-Ys(ii))^2);
        rho_dot = ((X-Xs(ii))*(Xdot - Xs_dot(ii)) + (Y - Ys(ii))*(Ydot - Ys_dot(ii))) / (rho);
        phi = phi(ii);
        
        yhat = [rho; rho_dot; phi];
        return
    end

end
end