%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Chris Chamberlain
% Written:  07 Dec 2017
% Revised:  07 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  ASEN 5044 - Statistical Estimation for Dynamical Systems Final
%           Project. The Kalman Filter equations were obtained from Simon, 
%           page ## (eq ##)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   states - state struct
%               .x - state vector
%               .dx - state perturbations vector
% 
%           inputs - input struct
%               .u - input vector 
%               .unom - nominal input vector 
% 
%           meas - measurement struct
%               .y - measurement vector
%               .dy - measurement perturbations vector
% 
%           G - control effect matrix
% 
%           Omega - 
% 
%           P - predicted state covariance matrix, P(0)
% 
%           Q - process noise covariance matrix
% 
%           R - measurement noise covariance matrix
% 
%           n - number of states
% 
%           tf - final time (simulation duration), s
% 
%           dt - time step, s
% 
% Outputs:  dxhat - predicted state vectors perturbations
% 
%           sigma - positive 2sigma bounds for all states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dxhat,sigma] = LinearizedKF(states,inputs,meas,G,P,Q,Omega,R,n,tf,dt,truth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dxhat,sigma] = LinearizedKF(states,inputs,meas,G,Omega,P,Q,R,n,tf,dt,mu)
xhat = states.x;        dxhat = states.dx;
u = inputs.u;           unom = inputs.unom;
y = meas.y;             ynom = meas.ynom;

% (nxn) identity matrix
I = eye(n);

% 2sigma error bounds
sigma(:,1) = 2*sqrt(diag(P));

% input pertrubations
du(:,1) = u(:,1) - unom(:,1);

% measurement perturbations
dy = y(:,1) - ynom(:,1);

for kk = 1:tf/dt
  Anom = Anominal(xhat(:,kk+1),mu);
  F = I + dt*Anom;
  
  % C and H?
  [y(:,kk+1), H] = measure(xnom(:,kk+1), kk, dt); % NONLINEAR MEASUREMENT
  
  % time update step (-) superscript
  dxhat(:,kk+1) = F*dxhat(:,kk)+G*du;
  P = F*P*F' + Omega*Q*Omega';
  du = u(:,kk+1) - unom(:,kk+1);
  K = P*H'*inv(H*P*H'+R);                                   % Kalman gain

  % measurement update step (+) superscript
  dxhat(:,kk+1) = dxhat(:,kk+1)+K*(dy-H*dxhat(:,kk+1)); % a posteriori
  P = (I-K*H)*P*(I-K*H)'+K*R*K';
  
  
  
  dy = y(:,kk+1) - ynom(:,kk+1);
  
  % 2sigma error bounds
  sigma(:,kk+1) = 2*sqrt(diag(P));
end

end
