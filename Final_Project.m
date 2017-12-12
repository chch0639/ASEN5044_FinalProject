%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chris Chamberlain and Mitchell Smith
% Written: 07 Dec 2017
% Revised: 12 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: ASEN 5044 - Statistical Estimation for Dynamical Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; plotsettings(14,2); close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
problem = 3;
plot_flag = 1;
save_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants
n = 4;                          % number of states
m = 2;                          % number of inputs
p = 3;                          % number of measurements
stations = 6;                   % number of tracking stations
mu = 398600;                    % km^3/s^2
omegaE = 2*pi/86400;            % rad/s
RE = 6378;                      % km
r0 = 6678;                      % km
dt = 10;                        % s
rng(100);
x0 = r0;                        % km
xdot0 = 0;                      % km/s
y0 = 0;                         % km
ydot0 = r0*sqrt(mu/r0^3);       % km/s
xnom = [x0,xdot0,y0,ydot0];     % initial state
x_init = xnom;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% data:
load Truth.mat
% Includes:
%   Qtrue       (2x2 covariance matrix)
%   Rtrue       (2x2 covariance matrix)
%   tvec        (vector of times)
%   ydata       (matrix of true measurements)
%   measLabels  (strings describing data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO
% the linear model is hitting check2 two timesteps before ode45. Wtf?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 1  % HW8 question 2
        %% part a
        % CT LTI
        Anom = Anominal(xnom,mu);
        
        % B matrix linearized around nominal
        Bnom = [0,0; 1,0; 0,0; 0,1];
        
        % C matrix linearized around nominal
        Cnom = zeros(p,n,stations);
        for ii = 1:stations
            % initial conditions for current station
            xs = [RE*cos((ii-1)*(pi/3)),0,RE*sin((ii-1)*(pi/3)),0];
            
            Cnom(:,:,ii) = Cnominal(xnom,xs);
        end
        clear xs ii
        
        %% part b
        % convert to DT LTI
        dt = 10;                    % s, given
        P = 2*pi*sqrt(r0^3/mu);
%         time = 0:dt:P+dt;
        time = 0:dt:14000;
        u = [0;0];
        deltax(:,1) = [0,0.075,0,-0.021];    % random values entered
        
        % calculate nominal trajectory
        [~,xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
        xnom = xnom';
        
        Fnom = eye(n) + dt*Anom;
        Gnom = dt*Bnom;
        
        Omeganom = [0,0; 1,0; 0,0; 0,1];
        
        %% part c
        % satellite visibilty checks
        check1 = @(theta,phi)theta-pi/2 <= phi && phi <= theta+pi/2;
        check2 = @(theta,phi)theta-pi/2 >= phi && phi <= theta+pi/2-2*pi;
        check3 = @(theta,phi)theta+2*pi-pi/2 <= phi && phi <= theta+pi/2+2*pi;
        
        deltax = zeros(n,length(time));
        deltax(:,1) = [0,0.075,0,-0.021];    % random values entered
        
        % initial conditions for each tracking station
        for ii = 1:stations
            theta0(ii) = (ii-1)*pi/3;
        end
        
        for kk = 1:length(time)
            % linearize A around current state
            Anom = Anominal(xnom(:,kk),mu);
            
            % Convert A to discrete time
            Fnom = eye(n) + dt*Anom;
            
            % state perturbations
            deltax(:,kk+1) = Fnom*deltax(:,kk) + Gnom*u;

            [deltay(:,kk), ~] = measure(xnom(:,kk), kk, dt, 'linear', deltax(:,kk+1));
            [y_nom(:,kk), ~] = measure(xnom(:,kk), kk, dt, 'nonlinear');
            y_linear(:,kk) = deltay(:,kk) + y_nom(:,kk);
            
%           % Uncomment these lines to use a nonlinear measurement for the linearized system!            
%             [y_linear(:,kk), ~] = measure(xnom(:,kk)+deltax(:,kk), kk, dt, 'nonlinear');
%             [y_nom(:,kk), ~] = measure(xnom(:,kk), kk, dt, 'nonlinear');
        end
        
        plot_time = 0:dt:(kk*dt);
        
        % DT nonlinear model
        x0 = x_init;
        deltax(:,1) = [0,0.075,0,-0.021];
        [TOUT,XOUT] = ode45(@(t,x)NLode(t,x,u,mu),time,x_init+deltax(:,1)',options);
        % [TOUT,XOUT] = ode45(@(t,x)NLode(t,x,u,mu),time,x_init,options);
        for kk = 1:length(TOUT)
            [y_nonlinear(:,kk), ~] = measure(XOUT(kk,:), kk, dt, 'nonlinear');
        end
        
        if plot_flag == 1
            y_str = {'$\delta_x$, km','$\delta_{\dot{x}}$, km/s','$\delta_y$, km',...
                '$\delta_{\dot{y}}$, km/s'};
            figure
            hold on; box on; grid on;
            suptitle('Part 1 -- State Perturbations')
            for ii = 1:n
                subplot(n,1,ii)
                hold on; box on; grid on;
                ylabel(y_str{ii})
                plot(plot_time, deltax(ii,:),'r')
                plot(TOUT', XOUT(:,ii)' - xnom(ii,:),'--b')
                % plot(TOUT', XOUT(:,ii)' - (xnom(ii,:) + deltax(ii,1:end-1)), 'r') - this was wrong
                if ii == 1
                    legend('Linearized', 'ODE45')
                end
            end
            if save_flag == 1
                drawnow
                printFigureToPdf('1StateErr', [8,8],'in');
            end
            
            
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km',...
                '$\dot{y}$, km/s'};
            figure
            hold on; box on; grid on;
            suptitle('States vs Time, Non-linear and Linearized Approximate Dynamics Simulation')
            for ii = 1:n
                subplot(n,1,ii)
                hold on; box on; grid on;
                ylabel(y_str{ii})
                plot(plot_time(1:end-1), deltax(ii,1:end-1)+xnom(ii,:),'r')
                plot(TOUT', XOUT(:,ii)','--b')
%                 xlim([plot_time(1) plot_time(end)])
%                 xlim([1600 2300])
                if ii == 1
                    legend('Linearized', 'ODE45')
                end
            end
            
            figure
            suptitle('Part 1 -- Measurements Over Time')
            subplot(3,1,1)
            hold on; box on; grid on;
            plot(time, y_linear(1,:), 'r')
            plot(TOUT', y_nonlinear(1,:), 'b--')
            legend('Linearized', 'ODE45')
            ylabel('$\rho$, km')
            subplot(3,1,2)
            hold on; box on; grid on;
            plot(time, y_linear(2,:), 'r')
            plot(TOUT', y_nonlinear(2,:), 'b--')
            ylabel('$\dot{\rho}$, km/s')
            subplot(3,1,3)
            hold on; box on; grid on;
            plot(time, y_linear(3,:), 'r')
            plot(TOUT', y_nonlinear(3,:), 'b--')
            ylabel('$\phi$, rad')
            xlabel('Time, s')


%             if save_flag == 1
%                 drawnow
%                 printFigureToPdf('1Meas', [8,8],'in');
%             end
%             
%             figure
%             suptitle('Part 1 -- Measurements Errors Over Time')
%             for ii = 1:stations
%                 subplot(3,1,1)
%                 hold on; box on; grid on;
%                 plot(time, deltay(3*ii-2,:), 'r')
%                 plot(TOUT', y(3*ii-2,:) - ynom(3*ii-2,:), 'b--')
%                 legend('Linearized', 'ODE45','Location','SouthWest')
%                 ylabel('$e_{\rho}$, km')
% %                 xlim([0 time(end)])
% %                 xlim([1600 2300])
%                 subplot(3,1,2)
%                 hold on; box on; grid on;
%                 plot(time, deltay(3*ii-1,:), 'r')
%                 plot(TOUT', y(3*ii-1,:) - ynom(3*ii-1,:), 'b--')
%                 ylabel('$e_{\dot{\rho}}$, km/s')
% %                 xlim([0 time(end)])
% %                 xlim([1600 2300])
%                 subplot(3,1,3)
%                 hold on; box on; grid on;
%                 plot(time, deltay(3*ii,:), 'r')
%                 plot(TOUT', y(3*ii,:) - ynom(3*ii,:), 'b--')
%                 ylabel('$e_{\phi}$, rad')
%                 xlabel('Time, s')
% %                 ylim([-0.06 0.06])
% %                 xlim([0 time(end)])
% %                 xlim([1600 2300])
%             end
%             if save_flag == 1
%                 drawnow
%                 printFigureToPdf('1MeasErr', [8,8],'in');
%             end
        end
        
    case 2  % linearized KF (LKF)
        tf = tvec(end);
        
        % initialize matrices
        states.x = zeros(n,length(tvec));
        states.x(:,1) = [r0; 0; 0; r0*sqrt(mu/r0^3)];
        states.dx = zeros(n,length(tvec));
        states.dx(:,1) = ones(4,1)*0.001;       % initial state perturbations
        inputs.u = zeros(m,length(tvec));
        inputs.unom = zeros(m,length(tvec));
        meas.y = zeros(p,length(tvec));
        meas.ynom = zeros(p,length(tvec));
        B = [0,0; 1,0; 0,0; 0,1];
        G = dt*B;
        Omega = [0,0; 1,0; 0,0; 0,1];
        u(:,1) = [0;0];
        P = 1e4*eye(n);
%         Q = 1e4*diag([0.1 0.01 0.1 0.01]);
%         R = diag([0.001 0.0001 0.1]);
        Q = Qtrue;
        R = Rtrue;
        
        % calculate nominal trajectory
        [tnom,states.xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
        
        [dxhat,sigma] = LinearizedKF(states,inputs,ydata,G,Omega,P,Q,R,n,tf,dt,mu);
        
        if plot_flag
            figure()
            hold on; box on; grid on; axis equal
            title('Satellite Orbit')
            plot(states.xnom(:,1)', states.xnom(:,3)', 'r')
            plot(dxhat(1,:) + states.xnom(:,1)', dxhat(3,:) +  states.xnom(:,3)', 'b')
            legend('Nominal', 'LKF')
            xlabel('X position')
            ylabel('Y position')
            
            figure()
            suptitle('States Over Time')
            subplot(2,1,1)
            hold on; grid on; box on;
            ylabel('X position')
            plot(tvec, dxhat(1,:) + states.xnom(:,1)', 'r')
            plot(tnom, states.xnom(:,1)', 'b')
            legend('LKF', 'Nominal', 'Location', 'Best')
            
            subplot(2,1,2)
            hold on; grid on; box on;
            ylabel('Y position')
            xlabel('Time [s]')
            plot(tvec, dxhat(3,:) + states.xnom(:,3)', 'r')
            plot(tnom, states.xnom(:,3)', 'b')
        end
        
        
    case 3  % extended KF (EKF)
        close all;
        
        tf = tvec(end);
        
        % State Vector Definition:
        %   x1: X   (x position)
        %   x2: X'  (x velocity)
        %   x3: Y   (y position)
        %   x4: Y'  (y velocity)
        
        % initial state:
        xhat0 = [r0; 0; 0; r0*sqrt(mu/r0^3)];
        
        % Inputs to EKF:
        u = [0;0];
        P0 = eye(4);
        Omega = [0,0; 1,0; 0,0; 0,1];
        
        %%%%%% Guesses for R and Q to simulate actual measurements %%%%%%%%
        Q = Qtrue;
        R = Rtrue;
        
        [tnom, xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);

        % Extended Kalman Filter
        [xhat,sigma, yhat] = EKF(xhat0,u,ydata,P0,Q,Omega,R,n,tf,dt);
        
        if plot_flag
            figure()
            hold on; grid on; box on; axis equal;
            title('Extended Kalman Filter')
            plot(xnom(:,1), xnom(:,3)', 'b')
            plot(xhat(1,:), xhat(3,:), 'r--')
            legend('Nominal', 'Filtered')
            
            figure()
            suptitle('States Over Time')
            subplot(2,1,1)
            hold on; box on; grid on;
            plot(tvec, xhat(1,:), 'r')
            plot(tnom, xnom(:,1), 'b--')
            ylabel('X position')
            legend('EKF', 'Nominal')
            subplot(2,1,2)
            hold on; box on; grid on;
            plot(tvec, xhat(3,:), 'r')
            plot(tnom, xnom(:,3), 'b--')
            ylabel('Y position')
            xlabel('Time [s]')
            
            figure()
            suptitle('Measurements Over Time')
            subplot(3,1,1)
            hold on; box on; grid on;
            plot(tvec, yhat(1,:), 'r')
            plot(tvec, ydata(1,:), 'b--')
            ylabel('$\rho$, km')
            legend('EKF Measurements', 'True Measurements', 'Location', 'Best')
            subplot(3,1,2)
            hold on; box on; grid on;
            plot(tvec, yhat(2,:), 'r')
            plot(tvec, ydata(2,:), 'b--')
            ylabel('$\dot{\rho}$, km/s')
            subplot(3,1,3)
            hold on; box on; grid on;
            plot(tvec, yhat(3,:), 'r')
            plot(tvec, ydata(3,:), 'b--')
            ylabel('$\phi$, rad')
            xlabel('Time [s]')
        end
        
        
        
        
        
    case 4  % estimate state trajectory for LKF and EKF
        
        
    otherwise
        error('Invalid problem number!');
end





