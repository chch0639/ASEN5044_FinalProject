%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chris Chamberlain and Mitchell Smith
% Written: 07 Dec 2017
% Revised: 18 Dec 2017
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
x0 = r0;                        % km
xdot0 = 0;                      % km/s
y0 = 0;                         % km
ydot0 = r0*sqrt(mu/r0^3);       % km/s
xnom = [x0,xdot0,y0,ydot0];     % initial state
x_init = xnom;                  % initial state
options = odeset('RelTol',1e-12,'AbsTol',1e-12);    % ode tolerance
Nsims = 1;                      % number of simulations for NEES and NIS
rng(100)                        % fix random number generator seed

% data:
load Truth.mat
% Includes:
%   Qtrue       (2x2 covariance matrix)
%   Rtrue       (2x2 covariance matrix)
%   tvec        (vector of times)
%   ydata       (matrix of true measurements)
%   measLabels  (strings describing data)
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
        u = [0;0];
        
        % calculate nominal trajectory
        [~,xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
        xnom = xnom';

        Fnom = eye(n) + dt*Anom;
        Gnom = dt*Bnom;
        
        Omeganom = [0,0; 1,0; 0,0; 0,1];
        
        %% part c
        % initial perturbations
        deltax = zeros(n,length(tvec));
        deltax(:,1) = [0,0.075,0,-0.021];
        
        % initial conditions for each tracking station
        for ii = 1:stations
            theta0(ii) = (ii-1)*pi/3;
        end
        
        for kk = 1:length(tvec)
            % linearize A around current state
            Anom = Anominal(xnom(:,kk),mu);
            
            % Convert A to discrete time
            Fnom = eye(n) + dt*Anom;
            
            % state perturbations
            deltax(:,kk+1) = Fnom*deltax(:,kk) + Gnom*u;

            %% check frst H for dy -- starting off with wrong rhodot
            [deltay(:,kk), ~] = measure(xnom(:,kk), kk, dt, 'linear', deltax(:,kk+1));
            [y_nom(:,kk), ~] = measure(xnom(:,kk), kk, dt, 'nonlinear');
             y_linear(:,kk) = deltay(:,kk) + y_nom(:,kk);
             
%           % Uncomment these lines to use a nonlinear measurement for the linearized system!            
%             [y_linear(:,kk), ~] = measure(xnom(:,kk)+deltax(:,kk), kk, dt, 'nonlinear');
%             [y_nom(:,kk), ~] = measure(xnom(:,kk), kk, dt, 'nonlinear');
        end

        % plot linear measurements as well as error compared to data
        figure
        title('Linear Measurement Errors')
        hold on; box on; grid on;
        plot(tvec(2:end),y_nom(1,2:end),'b')
        plot(tvec(2:end),y_linear(1,2:end),'r')
        ylabel('$\rho$, km')
        xlabel('Time, s')
        xlim([tvec(1) tvec(end)])

        % DT nonlinear model
        [TOUT,XOUT] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init+deltax(:,1)',options);

        for kk = 1:length(TOUT)
            [y_nonlinear(:,kk), ~] = measure(XOUT(kk,:), kk, dt, 'nonlinear');
        end

        %% plot results
        % state perturbations
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
                plot(tvec, deltax(ii,1:end-1),'r')
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
            
            % states vs time
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km',...
                '$\dot{y}$, km/s'};
            figure
            hold on; box on; grid on;
            suptitle('States vs Time, Non-linear and Linearized Approximate Dynamics Simulation')
            for ii = 1:n
                subplot(n,1,ii)
                hold on; box on; grid on;
                ylabel(y_str{ii})
                plot(tvec, deltax(ii,1:end-1)+xnom(ii,:),'r')
                plot(TOUT', XOUT(:,ii)','--b')
                xlim([tvec(1) tvec(end)])
                if ii == 1
                    legend('Linearized', 'ODE45')
                end
            end
       
            % measurements vs time
            figure
            suptitle('Part 1 -- Measurements Over Time')
            subplot(3,1,1)
            hold on; box on; grid on;
            plot(tvec(2:end), y_linear(1,2:end)-deltay(1,2:end), 'r')
%             plot(TOUT', y_nonlinear(1,:), 'b--')
            legend('Linearized', 'ODE45')
            ylabel('$\rho$, km')
            subplot(3,1,2)
            hold on; box on; grid on;
            plot(tvec(2:end), y_linear(2,2:end)-deltay(2,2:end), 'r')
%             plot(TOUT', y_nonlinear(2,:), 'b--')
            ylabel('$\dot{\rho}$, km/s')
            subplot(3,1,3)
            hold on; box on; grid on;
            plot(tvec(2:end), y_linear(3,2:end)-deltay(3,2:end), 'r')
%             plot(TOUT', y_nonlinear(3,:), 'b--')
            ylabel('$\phi$, rad')
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('1Meas', [8,8],'in');
            end
        end
        
        
    case 2  % linearized KF (LKF)
        tf = tvec(end);
        
        % initialize matrices
        states.x = zeros(n,length(tvec));
        states.x(:,1) = x_init;
        
        states.dx = zeros(n,length(tvec));
        states.dx(:,1) = [0,0.075,0,-0.021];  % initial state perturbations
        
        inputs.u = zeros(m,length(tvec));
        inputs.unom = zeros(m,length(tvec));
        
        meas.y = zeros(p,length(tvec));
        meas.ynom = zeros(p,length(tvec));
        
        B = [0,0; 1,0; 0,0; 0,1];
        G = dt*B;
        Gamma = B;
        Omega = dt*B;
        
        u(:,1) = [0;0];
        P = 1e6*eye(n);
        Q = Qtrue;
        SvQ = chol(Q)';
        R = Rtrue;
        SvR = chol(R)';
        
        for kk = 1:Nsims
          % Ccreate noisy states and measurements
          xnoise = x_init';
          for ii = 1:length(tvec)-1
            wtilde = SvQ*randn([length(Q) 1]);        % process noise
            
            % integrate one time step to find next initial state
            [~, x_temp] = ode45(@(t,x)NLode(t,x,u,mu, wtilde),[0 dt],xnoise(:,end),options);         
            xnoise(:,ii+1) = x_temp(end,:)';
            
            % use noisy state to find noisy measurement
            ynoise(:,ii+1) = measure(xnoise(:,ii+1), ii, dt, 'nonlinear');
          end

          for ii = 1:length(tvec)-1
              if ynoise(4,ii+1) ~= 0
                  vtilde = SvR*randn([length(R) 1]);  % measurement noise
                  ynoise(1:3,ii+1) = ynoise(1:3,ii+1) + vtilde;
              end
          end
          states.xnoise = xnoise;
          
          % calculate nominal trajectory
          [tnom,states.xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
          states.xnom = states.xnom';
          
          % linearized KF to get state perturbations
          [dxhat,dyhat,ynom, sigma, NEES(kk,:), NIS(kk,:)] =...
            LinearizedKF(states,inputs,ynoise,G,Omega,P,Q,R,n,tf,dt,mu);
        end
        
        % Perform NEES and NIS Tests
        NEESbar = mean(NEES,1);
        alpha_NEES = 0.05;
        Nnx = Nsims*n;
        r1x = chi2inv(alpha_NEES/2, Nnx)./Nsims;
        r2x = chi2inv(1-alpha_NEES/2, Nnx)./Nsims;
        
        NISbar = mean(NIS, 1);
        alpha_NIS = 0.05;
        Nny = Nsims*p;
        r1y = chi2inv(alpha_NIS/2, Nny)./Nsims;
        r2y = chi2inv(1-alpha_NIS/2, Nny)./Nsims;
        
        
        %% plot results
        if plot_flag
            figure
            hold on; box on; grid on; axis equal
            title('Linearized Kalman Filter Orbit')
            plot(states.xnoise(1,:), states.xnoise(3,:),'--k')
            plot(dxhat(1,:)+states.xnom(1,:),dxhat(3,:)+states.xnom(3,:),'r')
            legend('Noise','Estimated','Location','EastOutside')
            xlabel('X Position, km')
            ylabel('Y Position, km')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_orbit', [8,8],'in');
            end
            
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
            y_limits = [-530,530;-2,2;-300,300;-1,1];       
            figure
            suptitle('Linearized KF -- State Residuals')
            hold on; grid on; box on;
            for ii = 1:n
              subplot(4,1,ii)
              hold on; grid on; box on;
              plot(tvec,dxhat(ii,:),'r')
              plot(tvec,sigma(ii,:),'--k')
              if ii == 1
                legend('Residuals','2$\sigma$ Bounds','Location','SouthEast')
              end
              plot(tvec,-sigma(ii,:),'--k')
              ylim([y_limits(ii,1) y_limits(ii,2)])
              ylabel(y_str{ii})
            end
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_residuals', [8,8],'in');
            end
            
            figure
            suptitle('Linearized KF -- States Over Time')
            subplot(2,1,1)
            hold on; grid on; box on;
            plot(tvec, states.xnoise(1,:), 'k','Linewidth',4)
            plot(tvec, dxhat(1,:) + states.xnom(1,:), 'r')
            ylabel('X Position, km')
            legend('Noise', 'Estimated', 'Location', 'Best')
            subplot(2,1,2)
            hold on; grid on; box on;
            plot(tvec, states.xnoise(3,:), 'k','Linewidth',4)
            plot(tvec, dxhat(3,:) + states.xnom(3,:), 'r')
            ylabel('Y Position, km')
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_states', [8,8],'in');
            end
            
            y_str = {'$\rho$, km','$\dot{\rho}$, km/s','$\phi$, rad','Station ID'};
            figure()
            suptitle('Measurements Over Time')
            for ii = 1:p
              subplot(p+1,1,ii)
              hold on; box on; grid on;
              plot(tvec, ynoise(ii,:), 'k','Linewidth',4)
              plot(tvec, dyhat(ii,:)+ynom(ii,:), 'r')
              if ii == 1
                legend('Noise', 'Estimated','Location','Best')
              end
              ylabel(y_str{ii})
            end
            subplot(p+1,1,p+1)
            hold on; box on; grid on;
            plot(tvec,ynom(p+1,:),'r')
            ylabel(y_str{p+1})
            xlabel('Time, s') 
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_measurements', [8,8],'in');
            end
 
            figure()
            hold on; box on; grid on;
            title('NEES Test')
            h1 = plot(NEESbar, 'ro');
            h2 = plot(r1x*ones(size(NEESbar)), 'r--');
            plot(r2x*ones(size(NEESbar)), 'r--')
            xlabel('Time Step [k]')
            ylabel('NEES statistic, $\bar{\epsilon}_x$')
            ylim([0 r2x*1.5])
            legend([h1 h2], 'NEES @ time k', 'Bounds', 'Location', 'Best')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_NEES', [8,4],'in');
            end
        end
        
        
    case 3  % extended KF (EKF)
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
        P0 = 10*eye(4);
        Omega = dt*[0,0; 1,0; 0,0; 0,1];
        
        %%%%%% Guesses for R and Q to simulate actual measurements %%%%%%%%
        Q = Qtrue;
        R = Rtrue;
        
        % Create Nominal Conditions
        [tnom, xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
        
        for kk = 1:Nsims
            % Create Noisy Measurements
            xnoise = x_init';
            for ii = 1:length(tvec)-1
                Sv = chol(Q)';
                q = randn([length(Q) 1]);
                wtilde = Sv*q;
                
                [~, x_temp] = ode45(@(t,x)NLode(t,x,u,mu, wtilde),[0 dt],xnoise(:,end),options);
                xnoise(:,ii+1) = x_temp(end,:)';
                
                Sv = chol(R)';
                r = randn([length(R) 1]);
                ynoise(:,ii+1) = measure(xnoise(:,ii+1), ii, dt, 'nonlinear');
                
                if ynoise(4,ii+1) ~= 0
                    ynoise(1:3,ii+1) = ynoise(1:3,ii+1) + Sv*r;
                end
            end
            
            % Extended Kalman Filter
            [xhat,sigma, yhat, NEES(kk,:), NIS(kk,:)] =...
              EKF(x_init',xnoise,u,ynoise,P0,100*Q,Omega,R,n,tf,dt);  
        end
        
        % Perform NEES and NIS Tests
        NEESbar = mean(NEES,1);
        alpha_NEES = 0.05;
        Nnx = Nsims*n;
        r1x = chi2inv(alpha_NEES/2, Nnx)./Nsims;
        r2x = chi2inv(1-alpha_NEES/2, Nnx)./Nsims;
        
        NISbar = mean(NIS, 1);
        alpha_NIS = 0.05;
        Nny = Nsims*p;
        r1y = chi2inv(alpha_NIS/2, Nny)./Nsims;
        r2y = chi2inv(1-alpha_NIS/2, Nny)./Nsims;
        
        if plot_flag
          figure()
          hold on; grid on; box on; axis equal;
          title('Extended Kalman Filter Orbit')
          plot(xnoise(1,:), xnoise(3,:), 'k','Linewidth',4)
          plot(xhat(1,:), xhat(3,:), 'r--')
          xlabel('X Position, km')
          ylabel('Y Position, km')
          legend('Noise', 'Estimated')
          if save_flag == 1
            drawnow
            printFigureToPdf('EKF_orbit', [8,8],'in');
          end
          
          y_str = {'X Position, km','X Velocity, km/s','Y Position, km','Y Velocity, km/s'};
          figure()
          suptitle('EKF -- State Residuals')
          for ii = 1:n
            subplot(n,1,ii)
            hold on; box on; grid on;
            plot(xnoise(ii,:) - xhat(ii,:),'r')
            plot(sigma(ii,:), 'k--')
            if ii == 1
              legend('Residuals','2$\sigma$ Bounds','Location','SouthEast')
            end
            plot(-sigma(ii,:), 'k--')
            xlim([tvec(1)/dt tvec(end)/dt])
            if ii == 2 || ii == 4
              ylim([-0.5 0.5])
            end
            ylabel(y_str{ii})
          end
          xlabel('Time step, k')
          if save_flag == 1
            drawnow
            printFigureToPdf('EKF_residuals', [8,8],'in');
          end

          figure()
          suptitle('EKF -- States Over Time')
          subplot(2,1,1)
          hold on; box on; grid on;
          plot(tvec, xnoise(1,:), 'k','Linewidth',4)
          plot(tvec, xhat(1,:), 'r')
          xlim([tvec(1) tvec(end)])
          ylabel('X position, km')
          legend('Noise','Estimated','Location','Best')
          subplot(2,1,2)
          hold on; box on; grid on;
          plot(tvec, xnoise(3,:), 'k','Linewidth',4)
          plot(tvec, xhat(3,:), 'r')
          xlim([tvec(1) tvec(end)])
          ylabel('Y position, km')
          xlabel('Time, s')
          if save_flag == 1
            drawnow
            printFigureToPdf('EKF_states', [8,8],'in');
          end
          
          y_str = {'$\rho$, km','$\dot{\rho}$, km/s','$\phi$, rad','Station ID'};
          figure()
          suptitle('EKF -- Measurements Over Time')
          for ii = 1:p
            subplot(p+1,1,ii)
            hold on; box on; grid on;
            plot(tvec, ynoise(ii,:), 'k','Linewidth',4)
            plot(tvec, yhat(ii,:), 'r')
            xlim([tvec(1) tvec(end)])
            ylabel(y_str{ii})
            if ii == 1
              legend('Noise','Estimated','Location', 'Best')
            end
          end
          subplot(p+1,1,p+1)
            hold on; box on; grid on;
            plot(tvec, yhat(p+1,:), '-r')
            xlim([tvec(1) tvec(end)])
            ylabel(y_str{p+1})
          xlabel('Time, s')
          if save_flag == 1
            drawnow
            printFigureToPdf('EKF_measurements', [8,8],'in');
          end
          
          figure()
          hold on; box on; grid on;
          title('EKF -- NEES Test')
          h1 = plot(NEESbar, 'ro');
          h2 = plot(r1x*ones(size(NEESbar)), 'k--');
          plot(r2x*ones(size(NEESbar)), 'k--')
          xlabel('Time Step, k')
          ylabel('NEES Statistic, $\bar{\epsilon}_x$')
          % ylim([-r1y 1.5*r2y])
          xlim([tvec(1)/dt tvec(end)/dt])
          legend([h1 h2], 'NEES @ Time k', 'Bounds', 'Location', 'Best')
          if save_flag == 1
            drawnow
            printFigureToPdf('EKF_NEES', [8,4],'in');
          end
          
          figure()
          hold on; box on; grid on;
          title('EKF -- NIS Test')
          h1 = plot(NISbar, 'ro');
          h2 = plot(r1y*ones(size(NISbar)), 'k--');
          plot(r2y*ones(size(NISbar)), 'k--')
          xlabel('Time Step, k')
          ylabel('NIS Statistic, $\bar{\epsilon}_y$')
          % ylim([-r1y 1.5*r2y])
          xlim([tvec(1)/dt tvec(end)/dt])
          legend([h1 h2], 'NIS @ Time k', 'Bounds', 'Location', 'Best')
          if save_flag == 1
            drawnow
            printFigureToPdf('EKF_NIS', [8,4],'in');
          end
        end
        
         
    case 4  % estimate state trajectory for LKF and EKF
        
        
    otherwise
        error('Invalid problem number!');
end
