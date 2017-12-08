%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Mitchell Smith
% Written:  07 Dec 2017
% Revised:  07 Dec 2017
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
% Outputs:  yhat:   Measurement given state and input
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yhat = measure(xhat, u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function yhat = measure(xhat, k, dt)

stations = 6;           % Number of tracking stations
RE = 6378;              % Radius of Earth [km]
omegaE = 2*pi/86400;    % Rotation rate of Earth [rad/s]


% Extract States to make equations shorter and easier to read
X = xhat(1);
Xdot = xhat(2);
Y = xhat(3);
Ydot = xhat(4);

% Get state of each of the stations
for ii = 1:stations
    theta0 = (ii-1)*pi/3; % Initial position of each station
    
    % X and Y position and velocity of each staion
    Xs(ii) = RE*cos(omegaE*(k-1)*dt + theta0);
    Xs_dot(ii) = -RE*omegaE*sin(omegaE*(k-1)*dt + theta0);
    Ys(ii) = RE*sin(omegaE*(k-1)*dt + theta0);
    Ys_dot(ii) = RE*omegaE*cos(omegaE*(k-1)*dt + theta0);
    
    theta(ii) = atan2(Ys(ii)/Xs(ii));
    phi(ii) = atan2((Y-Ys(ii)), (X-Xs(ii)));
    if phi(ii) % TODO: get phi requirement
        
    end
    
end

rho(ii) = sqrt((X-Xs(ii))^2 + (Y-Ys(ii))^2);
rho_dot = ((X-Xdot)*() + ()*()) / ()

yhat = [rho; rho_dot; phi];

end