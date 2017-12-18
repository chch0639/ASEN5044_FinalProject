%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Mitchell Smith
% Written:  07 Dec 2017
% Revised:  18 Dec 2017
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
% Outputs:  xhat:   predicted total state vectors
% 
%           sigma:  positive 2sigma bounds for all states
% 
%           yhat:
% 
%           NEES:
% 
%           NIS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xhat,sigma, yhat, NEES, NIS] = EKF(xhat0,xtrue,u,y,P0,Q,Omega,R,n,tf,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xhat,sigma, yhat, NEES, NIS] = EKF(xhat0,xtrue,u,y,P0,Q,Omega,R,n,tf,dt)

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
  Anom = Anominal(xhat(:,kk), mu);
  F = I + dt*Anom;
    
  % time update step (-) superscript
  
  % Use full nonlinear dynamics to estimate the state
  dx = NLode([],xhat(:,kk),u,mu);
  
  % time update
  xhat(:,kk+1) = xhat(:,kk) + dt*dx;
  P = F*P*F' + Omega*Q*Omega'; % P-

  % measurement update step (+) superscript
  [yhat(:,kk+1), H] = measure(xhat(:,kk+1), kk, dt, 'nonlinear'); % NONLINEAR MEASUREMENT
  
  K = P*H'*inv(H*P*H'+R);                                   % Kalman gain
  
  % Compute NIS Statistic
  NIS(kk) = (y(1:3,kk+1)-yhat(1:3,kk+1))'*inv(H*P*H'+R)*(y(1:3,kk+1)-yhat(1:3,kk+1));
  
  if yhat(4,kk+1) == 0 || y(4,kk+1) ==0
      K = zeros(size(K));
      NIS(kk) = NaN;
  end
  xhat(:,kk+1) = xhat(:,kk+1)+K*(y(1:3,kk+1)-yhat(1:3,kk+1)); % a posteriori
  
  P = (I-K*H)*P; % P+
  
  % 2sigma error bounds
  sigma(:,kk+1) = 2*sqrt(diag(P));
  
  % Compute NEES statistic
  NEES(kk) = (xtrue(:,kk+1) - xhat(:,kk+1))'*inv(P)*(xtrue(:,kk+1) - xhat(:,kk+1));
   
end

end
