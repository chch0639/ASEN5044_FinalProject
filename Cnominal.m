%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chris Chamberlain and Mitchell Smith
% Written: 07 Nov 2017
% Revised: 07 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  ASEN 5044 - Statistical Estimation for Dynamical Final Project.
%           Calculate the continuous time (CT) C matrix for linearized
%           dynamics. Evaluate the C matrix for the given nominal
%           conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   xnom - nominal states
%           xs - tracking station states
% 
% Outputs:  Cnom - linearized C matrix around nominal conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cnom = Cnominal(x,xs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cnom = Cnominal(x,xs)
% 2D range, km
r = norm([x(1)-xs(1);x(3)-xs(3)]);

% build row 1
dh1_dx1 = (x(1)-xs(1))/r;
dh1_dx2 = 0;
dh1_dx3 = (x(3)-xs(3))/r;
dh1_dx4 = 0;

% build row 2
dh2_dx1 = ((x(2)-xs(2))/r) - ((x(1)-xs(1))*...
          ((x(1)-xs(1))*(x(2)-xs(2))+(x(3)-xs(3))*(x(4)-xs(4))))/r^3;
dh2_dx2 = (x(1)-xs(1))/r;
dh2_dx3 = ((x(4)-xs(4))/r) - ((x(3)-xs(3))*...
          ((x(1)-xs(1))*(x(2)-xs(2))+(x(3)-xs(3))*(x(4)-xs(4))))/r^3;
dh2_dx4 = (x(3)-xs(3))/r;

% build row 3
dh3_dx1 = (xs(3)-x(3))/((xs(1)-x(1))^2+(xs(3)-x(3))^2);
dh3_dx2 = 0;
dh3_dx3 = (x(1)-xs(1))/((xs(1)-x(1))^2+(xs(3)-x(3))^2);
dh3_dx4 = 0;

% combine for output
Cnom = [dh1_dx1,dh1_dx2,dh1_dx3,dh1_dx4;
        dh2_dx1,dh2_dx2,dh2_dx3,dh2_dx4;
        dh3_dx1,dh3_dx2,dh3_dx3,dh3_dx4];
end
