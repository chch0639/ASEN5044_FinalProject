%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chris Chamberlain and Mitchell Smith
% Written: 07 Dec 2017
% Revised: 19 Dec 2017
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
Nsims = 25;                     % number of simulations for NEES and NIS
rng(100)                        % fix random number generator seed
T = 3*2*pi*sqrt(r0^3/mu);       % three orbital periods

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
        
        Omeganom = dt*[0,0; 1,0; 0,0; 0,1];
        
        %% part c
        % initial perturbations
        deltax = zeros(n,length(tvec));
        deltax(:,1) = [0,0.075,0,-0.021];
        
        % initial conditions for each tracking station
        for ii = 1:stations
            theta0(ii) = (ii-1)*pi/3;
        end
        
        for kk = 1:length(tvec)-1
            % linearize A around current state
            Anom = Anominal(xnom(:,kk),mu);
            
            % Convert A to discrete time
            Fnom = eye(n) + dt*Anom;
            
            % state perturbations
            deltax(:,kk+1) = Fnom*deltax(:,kk) + Gnom*u;
            
            % linear and nominal measurements
            [y_linear(:,kk+1), ~] = measure(xnom(:,kk+1)+deltax(:,kk+1), kk, dt, 'nonlinear');
            [y_nom(:,kk+1), ~] = measure(xnom(:,kk+1), kk, dt, 'nonlinear');
        end
        
        % DT nonlinear model
        [TOUT,XOUT] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init+deltax(:,1)',options);
        
        for kk = 1:length(TOUT)-1
            [y_nonlinear(:,kk+1), ~] = measure(XOUT(kk+1,:), kk, dt, 'nonlinear');
        end
        
        %% plot results
        if plot_flag == 1
          % nominal and perturbed orbit comparison
          figure
          axis equal
          hold on; box on; grid on;
          title('Orbit Comparison')
          plot(xnom(1,:),xnom(3,:),'--k')
          plot(XOUT(:,1),XOUT(:,3),'r')
          legend('Nominal','Perturbed')
          xlabel('X Position, km')
          ylabel('Y Position, km')
          xlim([-6800 6800])
          ylim([-6800 6800])
          if save_flag == 1
            drawnow
            printFigureToPdf('Part1_orbit', [8,8],'in');
          end
          
          % state perturbations
            y_str = {'$\delta_x$, km','$\delta_{\dot{x}}$, km/s','$\delta_y$, km',...
                '$\delta_{\dot{y}}$, km/s'};
            figure
            hold on; box on; grid on;
            suptitle('Part 1 -- State Perturbations')
            for ii = 1:n
                subplot(n,1,ii)
                hold on; box on; grid on;
                ylabel(y_str{ii})
                plot(TOUT', XOUT(:,ii)' - xnom(ii,:),'k','Linewidth',4)
                plot(tvec, deltax(ii,:),'r')
                if ii == 1
                    legend('ode45','Linearized','Location','Best')
                end
            end
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('Part1_residuals', [8,8],'in');
            end
            
            
            % states vs time
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km',...
                '$\dot{y}$, km/s'};
            figure
            hold on; box on; grid on;
            suptitle('States Over Time')
            for ii = 1:n
                subplot(n,1,ii)
                hold on; box on; grid on;
                ylabel(y_str{ii})
                plot(TOUT', XOUT(:,ii)','k','Linewidth',4)
                plot(tvec, deltax(ii,:)+xnom(ii,:),'r')
                xlim([tvec(1) tvec(end)])
                if ii == 1
                    legend('ode45','Linearized','Location','Best')
                end
            end
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('Part1_states', [8,8],'in');
            end
            
            % measurements vs time
            figure
            suptitle('Part 1 -- Measurements Over Time')
            subplot(3,1,1)
            hold on; box on; grid on;
            plot(TOUT', y_nonlinear(1,:),'k','Linewidth',4)
            plot(tvec(2:end), y_linear(1,2:end), 'r')
            legend('ode45','Linearized','Location','Best')
            ylabel('$\rho$, km')
            subplot(3,1,2)
            hold on; box on; grid on;
            plot(TOUT', y_nonlinear(2,:),'k','Linewidth',4)
            plot(tvec(2:end), y_linear(2,2:end), 'r')
            ylabel('$\dot{\rho}$, km/s')
            subplot(3,1,3)
            hold on; box on; grid on;
            plot(TOUT', y_nonlinear(3,:),'k','Linewidth',4)
            plot(tvec(2:end), y_linear(3,2:end), 'r')
            ylabel('$\phi$, rad')
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('Part1_measurements', [8,8],'in');
            end
        end
        
        
    case 2  % linearized KF (LKF)
      tvec = 0:dt:T;  
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
        P0 = 100*diag([.01 .001 .01 .001]);
        Q = Qtrue;
        Q_LKF = Q;
        SvQ = chol(Q)';
        R = Rtrue;
        R_LKF = R;
        SvR = chol(R)';
        
        for kk = 1:Nsims
            % Create noisy states and measurements
            xnoise = x_init' + states.dx(:,1);
            for ii = 1:length(tvec)-1
                wtilde = SvQ*randn([length(Q) 1]);        % process noise
                
                % integrate one time step to find next initial state
                [~, x_temp] = ode45(@(t,x)NLode(t,x,u,mu, wtilde),[0 dt],xnoise(:,end),options);
                %[~, x_temp] = ode45(@(t,x)NLode(t,x,u,mu),[0 dt],xnoise(:,end),options);
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
            [dxhat,dyhat,ynom,ynoise, sigma,ex, NEES(kk,:), NIS(kk,:)] =...
                LinearizedKF(states,inputs,ynoise,G,Omega,P0,Q_LKF,R_LKF,n,tf,dt,mu);
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
            legend('Noise','Estimated')
            xlabel('X Position, km')
            ylabel('Y Position, km')
            xlim([-7500 7500])
            ylim([-7500 7500])
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_orbit', [8,8],'in');
            end
            
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
            y_limits = [-60,60;-0.25,0.25;-60,60;-0.25,0.25];
            figure
            suptitle('Linearized KF -- State Residuals')
            hold on; grid on; box on;
            for ii = 1:n
                subplot(4,1,ii)
                hold on; grid on; box on;
                plot(tvec,ex(ii,:),'r')
                plot(tvec,sigma(ii,:),'--k')
                if ii == 1
                    legend('Residuals','2$\sigma$ Bounds','Location','SouthEast')
                end
                plot(tvec,-sigma(ii,:),'--k')
                ylim([y_limits(ii,1) y_limits(ii,2)])
                xlim([tvec(1) tvec(end)])
                ylabel(y_str{ii})
            end
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_residuals', [8,8],'in');
            end
            
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
            figure
            suptitle('Linearized KF -- States Over Time')
            for ii = 1:n
              subplot(n,1,ii)
              hold on; grid on; box on;
              plot(tvec, states.xnoise(ii,:), 'k','Linewidth',4)
              plot(tvec, dxhat(ii,:) + states.xnom(ii,:), 'r')
              xlim([tvec(1) tvec(end)])
              ylabel(y_str{ii})
              if ii == 1 
                legend('Noise', 'Estimated', 'Location', 'Best')
              end
            end
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_states', [8,8],'in');
            end        
            
            y_str = {'$\rho$, km','$\dot{\rho}$, km/s','$\phi$, rad','Station ID'};
            figure()
            suptitle('Linearized KF -- Measurements Over Time')
            for ii = 1:p
                subplot(p+1,1,ii)
                hold on; box on; grid on;
                plot(tvec, ynoise(ii,:), 'k','Linewidth',4)
                plot(tvec, dyhat(ii,:)+ynom(ii,:), 'r')
                xlim([tvec(1) tvec(end)])
                if ii == 1
                    legend('Noise', 'Estimated','Location','Best')
                end
                ylabel(y_str{ii})
            end
            subplot(p+1,1,p+1)
            hold on; box on; grid on;
            plot(tvec,ynoise(p+1,:), 'k','Linewidth',4)
            plot(tvec,ynom(p+1,:),'r')
            xlim([tvec(1) tvec(end)])
            ylabel(y_str{p+1})
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_measurements', [8,8],'in');
            end
            
            figure()
            hold on; box on; grid on;
            title_string = sprintf('Linearized KF -- NEES Test: %d Trials', Nsims);
            title(title_string)
            h1 = plot(NEESbar, 'ro');
            h2 = plot(r1x*ones(size(NEESbar)), 'r--');
            plot(r2x*ones(size(NEESbar)), 'r--')
            xlabel('Time Step, k')
            ylabel('NEES statistic, $\bar{\epsilon}_x$')
%             ylim([0 15])
%             xlim([0 100])
            xlim([tvec(1)/dt tvec(end)/dt])
            legend([h1 h2], 'NEES @ time k', 'Bounds', 'Location', 'Best')
            if save_flag == 1
                drawnow
                printFigureToPdf('LKF_NEES', [8,4],'in');
            end
        end
        
        
    case 3  % extended KF (EKF)
      tvec = 0:dt:T;  
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
        P0 = diag([.01 .001 .01 .001]);

        Omega = dt*[0,0; 1,0; 0,0; 0,1];
        
        %%%%%% Guesses for R and Q to simulate actual measurements %%%%%%%%
        Q = Qtrue;
        Q_EKF = 0.95*Qtrue;
        R = Rtrue;
        R_EKF = Rtrue;
        
        % Create Nominal Conditions
        [tnom, xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
        
        for kk = 1:Nsims
            % Create Noisy Measurements
            x_vary = chol(P0)'*randn([length(P0) 1]);
            xnoise = x_init' + x_vary;
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
                EKF(x_init',xnoise,u,ynoise,P0,Q_EKF,Omega,R_EKF,n,tf,dt);
        end
        
        % Perform NEES and NIS Tests
        NEESbar = mean(NEES,1);
        alpha_NEES = 0.05;
        Nnx = Nsims*n;
        r1x = chi2inv(alpha_NEES/2, Nnx)./Nsims;
        r2x = chi2inv(1-alpha_NEES/2, Nnx)./Nsims;
        
        NISbar = nanmean(NIS, 1);
        alpha_NIS = 0.05;
        Nny = Nsims*p;
        r1y = chi2inv(alpha_NIS/2, Nny)./Nsims;
        r2y = chi2inv(1-alpha_NIS/2, Nny)./Nsims;
        
        %%
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
            
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
            y_limits = [-10 10;-0.075 0.075;-12 12;-0.075 0.075];
            figure()
            suptitle('EKF -- State Residuals')
            for ii = 1:n
                subplot(n,1,ii)
                hold on; box on; grid on;
                plot(xnoise(ii,:) - xhat(ii,:),'r')
                plot(sigma(ii,:), 'k--')
                if ii == 2
                    legend('Residuals','2$\sigma$ Bounds','Location','Best')
                end
                plot(-sigma(ii,:), 'k--')
                xlim([tvec(1)/dt tvec(end)/dt])
                ylim([y_limits(ii,1) y_limits(ii,2)])
                ylabel(y_str{ii})
            end
            xlabel('Time step, k')
            if save_flag == 1
                drawnow
                printFigureToPdf('EKF_residuals', [8,8],'in');
            end
            
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
            figure()
            suptitle('EKF -- States Over Time')
            for ii = 1:n
              subplot(n,1,ii)
              hold on; box on; grid on;
              plot(tvec, xnoise(ii,:), 'k','Linewidth',4)
              plot(tvec, xhat(ii,:), 'r')
              xlim([tvec(1) tvec(end)])
              ylabel(y_str{ii})
              if ii == 1
                legend('Noise','Estimated','Location','Best')
              end
            end
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
            plot(tvec, ynoise(p+1,:),'k','Linewidth',4)
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
            title_string = sprintf('EKF -- NEES Test: %d Trials', Nsims);
            title(title_string)
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
            title_string = sprintf('EKF -- NIS Test: %d Trials', Nsims);
            title(title_string)
            h1 = plot(NISbar, 'ro');
            h2 = plot(r1y*ones(size(NISbar)), 'k--');
            plot(r2y*ones(size(NISbar)), 'k--')
            xlabel('Time Step, k')
            ylabel('NIS Statistic, $\bar{\epsilon}_y$')
            ylim([0 10])
            xlim([tvec(1)/dt tvec(end)/dt])
            legend([h1 h2], 'NIS @ Time k', 'Bounds', 'Location', 'Best')
            if save_flag == 1
              drawnow
              printFigureToPdf('EKF_NIS', [8,4],'in');
            end
        end
        
        
    case 4  % estimate state trajectory for LKF and EKF
        %% Linearized KF
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
        P0 = 100*diag([.01 .001 .01 .001]);
        Q = Qtrue;
        Q_LKF = Q;
        SvQ = chol(Q)';
        R = Rtrue;
        R_LKF = R;
        SvR = chol(R)';
        
        % Create noisy states and measurements
        xnoise = x_init' + states.dx(:,1);
        for ii = 1:length(tvec)-1
          wtilde = SvQ*randn([length(Q) 1]);        % process noise
          
          % integrate one time step to find next initial state
          [~, x_temp] = ode45(@(t,x)NLode(t,x,u,mu, wtilde),[0 dt],xnoise(:,end),options);
          xnoise(:,ii+1) = x_temp(end,:)';
        end
        ynoise = ydata;
        states.xnoise = xnoise;
        
        % calculate nominal trajectory
        [~,states.xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
        states.xnom = states.xnom';
        
        % linearized KF to get state perturbations
        [dxhat,dyhat,ynom,ynoise, sigma, ex] =...
          LinearizedKF(states,inputs,ynoise,G,Omega,P0,Q_LKF,R_LKF,n,tf,dt,mu);
        
        
        %% plot results
        if plot_flag 
          y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
          figure
          suptitle('Linearized KF -- States Over Time')
          for ii = 1:n
            subplot(n,1,ii)
            hold on; grid on; box on;
            plot(tvec, dxhat(ii,:) + states.xnom(ii,:), 'r')
            plot(tvec, sigma(ii,:) + dxhat(ii,:) + states.xnom(ii,:),'--k')
            if ii == 1
              legend('States', '2$\sigma$ Bounds', 'Location', 'Best')
            end
            plot(tvec, -sigma(ii,:) + dxhat(ii,:) + states.xnom(ii,:),'--k')
            ylabel(y_str{ii})
          end
          xlabel('Time, s')
          if save_flag == 1
            drawnow
            printFigureToPdf('Part4_LKF_states', [8,8],'in');
          end
        end

        %% Extended KF
        clear xnoise ynoise xnom sigma xhat
        % State Vector Definition:
        %   x1: X   (x position)
        %   x2: X'  (x velocity)
        %   x3: Y   (y position)
        %   x4: Y'  (y velocity)
        
        % Inputs to EKF:
        u = [0;0];
        P0 = diag([.01 .001 .01 .001]);
        
        Omega = dt*[0,0; 1,0; 0,0; 0,1];
        
        %%%%%% Guesses for R and Q to simulate actual measurements %%%%%%%%
        Q = Qtrue;
        Q_EKF = 0.95*Qtrue;
        R = Rtrue;
        R_EKF = Rtrue;
        
        % Create Nominal Conditions
        [tnom, xnom] = ode45(@(t,x)NLode(t,x,u,mu),tvec,x_init,options);
        
        xnoise = x_init';
        for ii = 1:length(tvec)-1
          Sv = chol(Q)';
          q = randn([length(Q) 1]);
          wtilde = Sv*q;
          
          [~, x_temp] = ode45(@(t,x)NLode(t,x,u,mu, wtilde),[0 dt],xnoise(:,end),options);
          xnoise(:,ii+1) = x_temp(end,:)';
        end
        ynoise = ydata;
        
        % Extended Kalman Filter
        [xhat,sigma,yhat] = EKF(x_init',xnoise,u,ynoise,P0,Q_EKF,Omega,R_EKF,n,tf,dt);
        
        %%
        if plot_flag
          y_str = {'$\rho$, km','$\dot{\rho}$, km/s','$\phi$, rad','Station ID'};
          figure()
          suptitle('EKF -- Measurements Over Time')
          for ii = 1:n
            subplot(n,1,ii)
            hold on; box on; grid on;
            plot(tvec, ydata(ii,:), 'k','Linewidth',4)
            plot(tvec, yhat(ii,:), 'r')
            xlim([tvec(1) tvec(end)])
            ylabel(y_str{ii})
            if ii == 1
              legend('True','Estimated','Location','Best')
            end
          end
          xlabel('Time, s')
          if save_flag == 1
                drawnow
                printFigureToPdf('Part4_EKF_measurements', [8,8],'in');
            end
          
          y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
            figure()
            suptitle('EKF -- States Over Time')
            for ii = 1:n
              subplot(n,1,ii)
              hold on; box on; grid on;
              plot(tvec, xhat(ii,:), 'r')
              plot(tvec, sigma(ii,:)+xhat(ii,:), '--k')
              plot(tvec,-sigma(ii,:)+xhat(ii,:), '--k')
              if ii == 1
                legend('States','2$\sigma$ Bounds','Location','Best')
              end
              xlim([tvec(1) tvec(end)])
              ylabel(y_str{ii})
            end
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('Part4_EKF_states', [8,8],'in');
            end
            
            y_str = {'$x$, km','$\dot{x}$, km/s','$y$, km','$\dot{y}$, km/s'};
            figure()
            suptitle('States Comparison Over Time')
            for ii = 1:n
              subplot(n,1,ii)
              hold on; box on; grid on;
              plot(tvec, dxhat(ii,:) + states.xnom(ii,:), 'k','Linewidth',4)
              plot(tvec, xhat(ii,:), 'r')
              if ii == 1
                legend('LKF','EKF','Location','Best')
              end
              xlim([tvec(1) tvec(end)])
              ylabel(y_str{ii})
            end
            xlabel('Time, s')
            if save_flag == 1
                drawnow
                printFigureToPdf('Part4_EKF_statescomp', [8,8],'in');
            end
            
            
            
        end
        
    otherwise
        error('Invalid problem number!');
end
