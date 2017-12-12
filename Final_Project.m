%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chris Chamberlain and Mitchell Smith
% Written: 07 Dec 2017
% Revised: 09 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: ASEN 5044 - Statistical Estimation for Dynamical Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; %plotsettings(14,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
problem = 1;
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
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 1  % HW8 question 2
        x0 = r0;                        % km
        xdot0 = 0;                      % km/s
        y0 = 0;                         % km
        ydot0 = r0*sqrt(mu/r0^3);       % km/s
        xnom = [x0,xdot0,y0,ydot0];     % initial state
        x_init = xnom;
        
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
        
        % calculate nominal trajectory
        [~,xnom] = ode45(@(t,x)NLode(t,x,u,mu),time,xnom,options);
        xnom = xnom';
        
        Fnom = eye(n) + dt*Anom;
        Gnom = dt*Bnom;
        
        Omeganom = [0,0; 1,0; 0,0; 0,1];
        
        %% part c
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
            
            Hnom = [];
            for ii = 1:stations
                % initial conditions for current station
                Xs(ii,kk) = RE*cos(omegaE*(kk-1)*dt + theta0(ii));
                Xs_dot(ii,kk) = -RE*omegaE*sin(omegaE*(kk-1)*dt + theta0(ii));
                Ys(ii,kk) = RE*sin(omegaE*(kk-1)*dt + theta0(ii));
                Ys_dot(ii,kk) = RE*omegaE*cos(omegaE*(kk-1)*dt + theta0(ii));
                xs = [Xs(ii,kk); Xs_dot(ii,kk); Ys(ii,kk); Ys_dot(ii,kk)];
                
                % linearize C around current state
                Cnom = Cnominal(xnom(:,kk),xs);
                
                phi(ii) = atan2((deltax(3,kk+1)+xnom(3,kk) - Ys(ii,kk)), deltax(1,kk+1)+xnom(1,kk) - Xs(ii,kk));
                theta(ii) = atan2(Ys(ii,kk), Xs(ii,kk));
                
                if (theta(ii)-pi/2 <= phi(ii) && phi(ii) <= theta(ii)+pi/2) ||...
                   (theta(ii)-pi/2 >= phi(ii) && phi(ii) <= theta(ii)+pi/2-2*pi)||...
                   (theta(ii)+2*pi-pi/2 <= phi(ii) && phi(ii) <= theta(ii)+pi/2+2*pi)
                    Hnom = [Hnom; Cnom];
                else
                    Hnom = [Hnom; NaN*ones(size(Cnom))];
                end
            end
            
            deltay(:,kk) = Hnom*deltax(:,kk+1);
        end
        
        plot_time = 0:dt:(kk*dt);
        
        % DT nonlinear model
        x0 = x_init;
        
        options = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [TOUT,XOUT] = ode45(@(t,x)NLode(t,x,u,mu),time,xnom(:,1)+deltax(:,1),options);
        
        y = NaN*ones(stations*3,length(TOUT));
        ynom = NaN*ones(stations*3,length(TOUT));
        for kk = 1:length(TOUT)
            for ii = 1:stations
                
                % initial conditions for current station
                Xs(ii,kk) = RE*cos(omegaE*(kk-1)*dt + theta0(ii));
                Xs_dot(ii,kk) = -RE*omegaE*sin(omegaE*(kk-1)*dt + theta0(ii));
                Ys(ii,kk) = RE*sin(omegaE*(kk-1)*dt + theta0(ii));
                Ys_dot(ii,kk) = RE*omegaE*cos(omegaE*(kk-1)*dt + theta0(ii));
                
                xs = [Xs(ii,kk); Xs_dot(ii,kk); Ys(ii,kk); Ys_dot(ii,kk)];
                
                phi = atan2((XOUT(kk,3) - Ys(ii,kk)), XOUT(kk,1) - Xs(ii,kk));
                theta = atan2(Ys(ii,kk), Xs(ii,kk));
                
                if (theta-pi/2 <= phi && phi <= theta+pi/2) ||...
                   (theta-pi/2 >= phi && phi <= theta+pi/2-2*pi)||...
                   (theta+2*pi-pi/2 <= phi && phi <= theta+pi/2+2*pi)
                      
                    X = XOUT(kk,1);
                    Xdot = XOUT(kk,2);
                    Y = XOUT(kk,3);
                    Ydot = XOUT(kk,4);
                    
                    rho = norm([X-xs(1); Y-xs(3)]);
                    rho_dot = ((X-xs(1))*(Xdot - xs(2)) + (Y-xs(3))*(Ydot-xs(4)))/(rho);
                    y(3*ii-2:3*ii,kk) = [rho; rho_dot; phi];
                end
                
                phi = atan2((xnom(3,kk) - Ys(ii,kk)), xnom(1,kk) - Xs(ii,kk));
                if (theta-pi/2 <= phi && phi <= theta+pi/2) ||...
                   (theta-pi/2 >= phi && phi <= theta+pi/2-2*pi)||...
                   (theta+2*pi-pi/2 <= phi && phi <= theta+pi/2+2*pi)
                    X = xnom(1,kk);
                    Xdot = xnom(2,kk);
                    Y = xnom(3,kk);
                    Ydot = xnom(4,kk);
                    
                    rho = norm([X-xs(1); Y-xs(3)]);
                    rho_dot = ((X-xs(1))*(Xdot - xs(2)) + (Y-xs(3))*(Ydot-xs(4)))/rho;
                    ynom(3*ii-2:3*ii,kk) = [rho; rho_dot; phi];
                end
            end
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
                xlim([plot_time(1) plot_time(end)])
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
                xlim([plot_time(1) plot_time(end)])
                if ii == 1
                    legend('Linearized', 'ODE45')
                end
            end
            
            
            
            figure
            suptitle('Part 1 -- Measurements Over Time')
            for ii = 1:stations
                subplot(3,1,1)
                hold on; box on; grid on;
                plot(time, deltay(3*ii-2,:) + ynom(3*ii-2,:), 'r')
                plot(TOUT', y(3*ii-2,:), 'b--')
                legend('Linearized', 'ODE45')
                ylabel('$\rho$, km')
                xlim([0 time(end)])
                subplot(3,1,2)
                hold on; box on; grid on;
                plot(time, deltay(3*ii-1,:) + ynom(3*ii-1,:), 'r')
                plot(TOUT', y(3*ii-1,:), 'b--')
                ylabel('$\dot{\rho}$, km/s')
                xlim([0 time(end)])
                subplot(3,1,3)
                hold on; box on; grid on;
                plot(time, deltay(3*ii,:) + ynom(3*ii,:), 'r')
                plot(TOUT', y(3*ii,:), 'b--')
                ylabel('$\phi$, rad')
                xlabel('Time, s')
                xlim([0 time(end)])
            end
            if save_flag == 1
                drawnow
                printFigureToPdf('1Meas', [8,8],'in');
            end
            
            figure
            suptitle('Part 1 -- Measurements Errors Over Time')
            for ii = 1:stations
                subplot(3,1,1)
                hold on; box on; grid on;
                plot(time, deltay(3*ii-2,:), 'r')
                plot(TOUT', y(3*ii-2,:) - ynom(3*ii-2,:), 'b--')
                legend('Linearized', 'ODE45','Location','SouthWest')
                ylabel('$e_{\rho}$, km')
                xlim([0 time(end)])
                subplot(3,1,2)
                hold on; box on; grid on;
                plot(time, deltay(3*ii-1,:), 'r')
                plot(TOUT', y(3*ii-1,:) - ynom(3*ii-1,:), 'b--')
                ylabel('$e_{\dot{\rho}}$, km/s')
                xlim([0 time(end)])
                subplot(3,1,3)
                hold on; box on; grid on;
                plot(time, deltay(3*ii,:), 'r')
                plot(TOUT', y(3*ii,:) - ynom(3*ii,:), 'b--')
                ylabel('$e_{\phi}$, rad')
                xlabel('Time, s')
%                 ylim([-0.06 0.06])
                xlim([0 time(end)])
            end
            if save_flag == 1
                drawnow
                printFigureToPdf('1MeasErr', [8,8],'in');
            end
        end
        
    case 2  % linearized KF (LKF)
        T = 2*pi*sqrt(r0^3/mu);           % period, s
        time = 0:dt:T+dt;                 % time vector, s
        tf = time(end);                   % final time, s
        
        % initialize matrices
        states.x = zeros(n,length(time));
        states.x(:,1) = [r0; 0; 0; r0*sqrt(mu/r0^3)];
        states.dx = zeros(n,length(time));
        states.dx(:,1) = ones(4,1)*0.001;       % initial state perturbations
        inputs.u = zeros(m,length(time));
        inputs.unom = zeros(m,length(time));
        meas.y = zeros(p,length(time));
        meas.ynom = zeros(p,length(time));
        B = [0,0; 1,0; 0,0; 0,1];
        G = dt*B;
        Omega = [0,0; 1,0; 0,0; 0,1];
        u(:,1) = [0;0];
        P = 1e4*eye(n);
        Q = 1e4*diag([0.1 0.01 0.1 0.01]);
        R = diag([0.001 0.0001 0.1]);
        Sv = chol(Q)';
        
        % calculate nominal trajectory
        states.x = [r0.*cos(2*pi.*time./T);
                    -2*pi*r0/T.*sin(2*pi.*time./T);
                    r0.*sin(2*pi.*time./T);
                    2*pi*r0/T.*cos(2*pi.*time./T)];
                  
        % CHECK THIS -- which is x, dx, y, ynom          
%         for kk = 1:length(time)
%             % Simulate noisy measurements and state perturbations
%             states.dx(:,kk) = Sv*randn(n,1);
%             meas.y(:,kk) = measure(states.dx(:,kk), kk, dt);
%         end          
        
        
        [dxhat,sigma] = LinearizedKF(states,inputs,meas,G,Omega,P,Q,R,n,tf,dt,mu);
        
        
    case 3  % extended KF (EKF)
        close all;
        
        T = 2*pi*sqrt(r0^3/mu);           % period, s
        time = 0:dt:T+dt;                 % time vector, s
        tf = time(end);                   % final time, s
        
        % State Vector Definition:
        %   x1: X   (x position)
        %   x2: X'  (x velocity)
        %   x3: Y   (y position)
        %   x4: Y'  (y velocity)
        
        % initial state:
        xhat0 = [r0; 0; 0; r0*sqrt(mu/r0^3)];
        
        % Inputs to EKF:
        u = [0;0];
        P0 = 1e4*eye(4);
        Omega = eye(4);
        
        %%%%%% Guesses for R and Q to simulate actual measurements %%%%%%%%
        Q = 1e4*diag([0.1 0.01 0.1 0.01]);
        Sv = chol(Q);
        R = diag([0.001 0.0001 0.1]);
        
        xnom = [r0.*cos(2*pi.*time./T);
            -2*pi*r0/T.*sin(2*pi.*time./T);
            r0.*sin(2*pi.*time./T);
            2*pi*r0/T.*cos(2*pi.*time./T)];
        
        y = zeros(3,length(time));
        for kk = 1:length(time)
            % Simulate noisy measurements
            q = randn([n 1]);
            xnom(:,kk) = xnom(:,kk) + Sv*q;
            y(:,kk) = measure(xnom, kk, dt);
        end
        
        % Extended Kalman Filter
        [xhat,sigma] = EKF(xhat0,u,y,P0,Q,Omega,R,n,tf,dt);
        
        figure()
        hold on; grid on; box on; axis equal;
        title('Extended Kalman Filter')
        plot(xnom(1,:), xnom(3,:), 'b')
        plot(xhat(1,:), xhat(3,:), 'r')
        legend('Noisy Data', 'Filtered')
        
        
        
        
        
        
    case 4  % estimate state trajectory for LKF and EKF
        
        
    otherwise
        error('Invalid problem number!');
end





