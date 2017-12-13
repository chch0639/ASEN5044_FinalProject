%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chris Chamberlain and Mitchell Smith
% Written: 07 Dec 2017
% Revised: 07 Dec 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:  ASEN 5044 - Statistical Estimation for Dynamical Final Project.
%           Calculate the continuous time (CT) A matrix for linearized
%           dynamics. Evaluate the A matrix for the given nominal
%           conditions. The dynamics within this matrix are based off of
%           given information for HW07.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   xnom - nominal states (nx1)
% 
%           mu - gravitational parameter
% 
% Outputs:  Anom - linearized A matrix around nominal conditions (nxn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anom = Anominal(xnom,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Anom = Anominal(xnom,mu)
% 2D range, km
r = norm([xnom(1);xnom(3)]);

% build row 1
df1_dx1 = 0;
df1_dx2 = 1;
df1_dx3 = 0;
df1_dx4 = 0;

% build row 2
df2_dx1 = mu*((2*xnom(1)^2-xnom(3)^2)/r^5);
df2_dx2 = 0;
df2_dx3 = 3*mu*xnom(1)*xnom(3)/r^5;
df2_dx4 = 0;

% build row 3
df3_dx1 = 0;
df3_dx2 = 0;
df3_dx3 = 0;
df3_dx4 = 1;

% build row 4
df4_dx1 = 3*mu*xnom(1)*xnom(3)/r^5;
df4_dx2 = 0;
df4_dx3 = mu*((2*xnom(3)^2-xnom(1)^2)/r^5);
df4_dx4 = 0;

% combine for output
Anom = [df1_dx1,df1_dx2,df1_dx3,df1_dx4;
        df2_dx1,df2_dx2,df2_dx3,df2_dx4;
        df3_dx1,df3_dx2,df3_dx3,df3_dx4;
        df4_dx1,df4_dx2,df4_dx3,df4_dx4];     
end
