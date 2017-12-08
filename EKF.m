%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Mitchell Smith
% Written:  07 Dec 2017
% Revised:  07 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  ASEN 5044 - Statistical Estimation for Dynamical Systems Final
%           Project. The Kalman Filter equations were obtained from Simon, 
%           page ## (eq ##)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   xhat0:  Initial estimate for state [nx1]
% 
%           u:          Input vector
%
%           y:          Measurement Vector
% 
%           Omega:  
% 
%           P0:         Predicted initial state covariance matrix, P(0)
% 
%           Q:          Process noise covariance matrix
% 
%           R:          Measurement noise covariance matrix
% 
%           n:          Number of states
% 
%           tf:         Final time (simulation duration), s
% 
%           dt:         Time step, s
%
%           Anominal:   Function for evaluating A matrix at the nominal
%                       state
%           
%           NL_eval:    Function for evaluating the nonlinear dynamics of
%                       the system. Takes the same form as an ODE function
%
%           NL_meas:    Function for evaluating the measurement of the
%                       dynamics of the system
% 
% Outputs:  xhat:   predicted total state vectors
% 
%           sigma - positive 2sigma bounds for all states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xhat,sigma] = LinearizedKF(xhat0,input,meas,P,Q,Omega,R,n,tf,dt, @Anominal, @NL_eval, @NL_meas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xhat,sigma] = EKF(xhat0,u,y,P0,Q,Omega,R,n,tf,dt, Anominal, NL_eval, NL_meas)

% Standard Gravitational Parameter
mu = 398600;

% Initial guesses become first state and covariance
xhat = xhat0;
P = P0;

% (nxn) identity matrix
I = eye(n);

% 2sigma error bounds
sigma(:,1) = 2*sqrt(diag(P));

for kk = 1:tf/dt
  Anom = Anominal(xhat(:,kk+1), mu);
  F = I + dt*Anom;
  
  % C and H?
  
  % time update step (-) superscript
  
  % Use full nonlinear dynamics to estimate the state used before
  % measurement update
  dx = NL_eval([],xhat(:,kk),u,mu);
  
  xhat(:,kk+1) = xhat(:,kk) + dt*dx;
  P = F*P*F' + Omega*Q*Omega';
  K = P*H'*inv(H*P*H'+R);                                   % Kalman gain

  % measurement update step (+) superscript
  yhat(:,kk+1) = NL_meas(xhat(:,kk+1), kk, dt); % NONLINEAR MEASUREMENT
  xhat(:,kk+1) = xhat(:,kk+1)+K*(y-yhat(:,kk+1)); % a posteriori
  P = (I-K*H)*P; % TODO: check this equation
  
  % 2sigma error bounds
  sigma(:,kk+1) = 2*sqrt(diag(P));
end

end
