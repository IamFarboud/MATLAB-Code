UTS = 73.98E6;             % Modulus of Rupture of Cypress (Pa)
%% 
b = 1.08;               % wing span (m)
l = 0.0066;             % width of the spar's cross-section (m)
A = l^2;                % Area of the spar's cross-section (m^2)
t_max = 0.02603;        % Max thickness of the Aerofoil (m)
y = t_max/2;            % distance from nutrual axis (m)
h = (t_max/2) - (l/2);  % distance from neutral axis and axis of spar's cross-section (m)
%% 
m = 1.8;                % Total mass (kg)
W = m*9.81;             % Total weight(N)
n = 2;                  % Load Factor
n_ult=n*1.5;            % Ultimate Load Factor
%% 
L = (W*n)/2;            % Lift force on each side (N)
L_ult = (W*n_ult)/2;    % Ultimate lift force on each side (N)
M_max = L*(b/2);        % Maximum bending moment at the tip chord (N/m)
M_max_ult = L_ult*b/4;  % Ultimate maximum bending moment at the tip chord (N/m)
%% 
I_gg = l^4/12;          % second moment of area of spar about its axis (m^4)
I_xx = I_gg + A*h^2;    % second moment of area of spar about neutral axis (m^4)
I_tot = 2*I_xx;         % second moment of area of both spars about neutral axis (m^4)
%% 
sigma_max = M_max*y/I_tot;          % Max bending stress on one side of the wing (Pa)
sigma_max_ult = M_max_ult*y/I_tot;  % Ultimate max bending stress on one side of the wing (Pa)
SF = UTS/sigma_max;                 % Safety Factor
SF_ult = UTS/sigma_max_ult;         % Ultimate safety Factor
