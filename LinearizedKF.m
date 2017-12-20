%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Chris Chamberlain and Mitchell Smith
% Written:  07 Dec 2017
% Revised:  18 Dec 2017
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
%           dy - 
% 
%           ynom - 
%
%           sigma - positive 2sigma bounds for all states
% 
%           NEES - 
% 
%           NIS - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dxhat,dy,ynom,sigma,NEES,NIS] = LinearizedKF(states,inputs,ydata,G,Omega,P,Q,R,n,tf,dt,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dxhat,dy,ynom,ydata,sigma,NEES,NIS] = LinearizedKF(states,inputs,ydata,G,Omega,P,Q,R,n,tf,dt,mu)
xnom = states.xnom;     dxhat = states.dx;      xtrue = states.xnoise;
u = inputs.u;           unom = inputs.unom;

% (nxn) identity matrix
I = eye(n);

% 2sigma error bounds
sigma(:,1) = 2*sqrt(diag(P));

% input pertrubations
du(:,1) = u(:,1) - unom(:,1);

for kk = 1:tf/dt
    % F_k
    Anom = Anominal(xnom(:,kk),mu);
    F = I + dt*Anom;
    
    % time update step (-) superscript
    dxhat(:,kk+1) = F*dxhat(:,kk)+G*du;
    P = F*P*F' + Omega*Q*Omega';
    
    have_measurement = ydata(4,kk+1) ~= 0;
    if have_measurement
        [ynom(:,kk+1), H] = measure(xnom(:,kk+1), kk, dt, 'nonlinear');
        have_measurement = have_measurement && ynom(4,kk+1) ~= 0;
        K = P*H'*inv(H*P*H'+R);
    else
        ynom(:,kk+1)= zeros(4,1);
    end
    
    if ~have_measurement
        ydata(:,kk+1) = zeros(4,1);
    end
    
    du = u(:,kk+1) - unom(:,kk+1);
    
    if have_measurement
        dy(:,kk+1) = ydata(1:3,kk+1) - ynom(1:3,kk+1);

        % Compute NIS Statistic
        NIS(kk) = dy(1:3,kk+1)'*inv(H*P*H'+R)*dy(1:3,kk+1);
    else
        dy(:,kk+1) = zeros(size(ydata(1:3, kk+1)));
        NIS(kk) = NaN;
    end
    
    if have_measurement
        % measurement update step (+) superscript
        dxhat(:,kk+1) = dxhat(:,kk+1)+K*(dy(1:3,kk+1)-H*dxhat(:,kk+1)); % a posteriori
        P = (I-K*H)*P;
    end
    
    % 2sigma error bounds
    sigma(:,kk+1) = 2*sqrt(diag(P));
    
    % Compute NEES statistic
    ex = xtrue(:,kk+1) - (xnom(:,kk+1) + dxhat(:,kk+1));
    NEES(kk) = ex'*inv(P)*ex;
end

end
