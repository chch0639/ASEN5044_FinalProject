%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chris Chamberlain and Mitchell Smith
% Written: 22 Nov 2017
% Revised: 05 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: ASEN 5044 - Statistical Estimation for Dynamical Systems HW 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx] = NLode(t,x,u,mu, varargin)

if isempty(varargin)
    wtilde = zeros([2 1]);
else
    wtilde = varargin{1};
end

r = sqrt(x(1)^2 + x(3)^2);

dx(1) = x(2);
dx(2) = -mu*x(1)/r^3 + u(1) + wtilde(1);
dx(3) = x(4);
dx(4) = -mu*x(3)/r^3 + u(2) + wtilde(2);

dx = dx';
end
