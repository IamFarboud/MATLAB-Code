clear
% Wing section airfoil name: WORTMANN FX 63-137 AIRFOIL (fx63137-il)   
b = 1.08;                     % wing span in meter
c= 0.19;    % wing chord
over_fuse = 0.15*c;  %fuselage overlap with wing
over_tp = 0.00539575; % fuselage overlap with tailplane
AR = b/c;                  % wing aspect ratio
S= b*c;                   % wing area
Snet = S-over_fuse;     % Wing area (net) (omit the overlap in fuselage)
c_tp = 0.12;    % chord of the tail plane
b_tp = 0.35;   % span of the tailplane
t_board = 0.005; % thickness of the foam board
rho_foam = 38.15; %density of foam kg/m^3
rho_board = 94.862; %density of foam board kg/m^3
CSA_wing = 0.00300456; % cross sectional area of the wing
CSA_tp = t_board*c_tp; % cross sectional area of the tailplane
CSA_fuse = 0.0225; % cross sectional area of the fuselage
CSA_undercut = 0.00157817; % cross sectional area of the wing undercut
rho_air = 1.225; % density of air kg/m^3
v_cruise = 13; % cruise velocity
mu_air = 0.000017894; % viscosity of air
Cl_max2D = 1.65; % max lift coefficient for infinite span wing
Cl0_2D = 0.8777; % zero angle of attack lift coefficient for infinite span wing
a2D = 0.0987*360/(2*pi); % lift slope curve for infinite span
AR_tp = b_tp/c_tp; % aspect ratio of tailplane
S_tp = b_tp*c_tp; % planform area of the tailplane
Snet_tp = S_tp - over_tp; % net area of the tailplane
St_S = Snet_tp/S; % proportion of the net tail area to the wing area

% --------------------------------------------------------
W_tennis_balls = 10 * 0.058*9.81;  % the weight of the tennis balls in N
W_speed_controller = 0.061*9.81;   % the weight of the speed controller in N
W_Servos = 4 * 0.014 * 9.81;       % the weight of the speed controller in N
W_battery = 0.175 *9.81;           % the weight of the battery in N
W_motor = 0.17*9.81;               % the weight of the motor in N
W_wheels = 3*0.023*9.81;           % the weight of the wheels in N
W_Nose_leg = 0.018*9.81;              % the weight of the battery in N
W_undercarriage = 0.12*9.81;       % the weight of the undercarriage  in N
W_wing = rho_foam*CSA_wing*b*9.81; % the weight of thew wing in N
W_tail = rho_board*t_board*c_tp*b_tp*9.81; % the weigh of the tail in N
W_fuse = 0.14025*9.81; % the wight of the fuselage in N
W_fin = rho_board*t_board*0.0078193*9.81; % the weight of the fin in N
W_undercut = rho_foam*CSA_undercut*0.15*9.81; % the weight of the undercut for the wing

% Weight and mass for max payload and no payload
W_max = W_tennis_balls + W_speed_controller + W_Servos + W_battery + ...
    + W_motor + W_wheels + W_Nose_leg + W_undercarriage + W_wing + W_tail + W_fuse + ...
    W_fin + W_undercut; 
mass_max=W_max/9.81; % total mass
W_NET = W_max - W_tennis_balls;
mass_NET = W_NET/9.81;
% --------------------------------------------------------
Re_cruise = (rho_air*v_cruise*c)/mu_air; % Reynolds number of wing in cruise
e = 1.78*(1-0.045*AR^0.68)-0.64; % oswalds coefficient for rectangular wings
a3D = a2D/(1+a2D/(pi*e*AR)); % lift curve slope for finite span wing
Cl_max3D = Cl_max2D*(a3D/a2D); % estimate of the maximum lift coefficient for finite span wing
% -----------------------------------------------------------
% Drag Calculations
Swet_wing = 0.4; % wetted area of the wing
Swet_fin = 0.01648706; % wetted area of the fin
Swet_tp = 0.0887; % wetted area of the tailplane
Swet_fuse = 0.31069336; % wetted area of the fuselage
Swet_wheel = 0.0420473; % wetted area of the fuselage
Swet_leg = 0.00070215; % wetted area of the front leg 
Swet_strut = 0.01995; % wetted area of the strut
Swet_sf = Swet_wing + Swet_fin + Swet_tp + Swet_fuse; % total wetted skin friction area
l_fin = c_tp; % characteristic length of the fin
l_fuse = 0.67; % characteristic length of the fuselage
Re_fuse = (rho_air*v_cruise*l_fuse)/mu_air; % reynolds number for skin friction of fuselage
Re_fin = (rho_air*v_cruise*l_fin)/mu_air; % reynolds number for skin friction of fin
Re_tp = (rho_air*v_cruise*c_tp)/mu_air; % reynolds number for skin friction of the tailplane
int = 1.3; % Assumed factor for the interference drag
cf_wing = 0.455/(log10(Re_cruise))^2.58; % skin friction coefficient of the wing
cf_fin = 0.455/(log10(Re_fin))^2.58; % skin friction coefficient of the fin
cf_tp = 0.455/(log10(Re_tp))^2.58; % skin friction coefficient of the tailplane
cf_fuse = 0.455/(log10(Re_fuse))^2.58; % skin friction coefficient of the fuselage
cf_total = 1.3*(cf_wing + cf_tp + cf_fuse + cf_fin); % total skin friction coefficient
cd_fuse = 2.1; % form drag coefficient of the fuselage
cd_wing = 0.4; % form drag coefficient of the wing
cd_tail = 2; % form drag coefficient of the tailplane
cd_wheel = 0.3; % form drag coefficient of the wheel
cd_strut = 2; % form drag coefficient of the strut
cd_leg = 0.3; % form drag coefficient for the front leg
Cd0 = int*(cf_fuse*(Swet_fuse/Swet_wing) + cf_wing + cf_fin*(Swet_fin/Swet_wing) + ...
    cf_tp*(Swet_tp/Swet_wing) + cd_fuse*(CSA_fuse/Swet_wing) + cd_wing*(CSA_wing/Swet_wing) ...
    + cd_tail*(CSA_tp/Swet_wing) + cd_wheel*3*(Swet_wheel/Swet_wing) + ...
    cd_strut*(Swet_strut/Swet_wing) + cd_leg*2*(Swet_leg/Swet_wing));

k = 1/(pi*AR*e); % induced drag constant
vs1 = sqrt(W_max/(0.5*rho_air*S*Cl_max3D)); % stall speed for max payload weight
vto1 = 1.2*vs1; % take off speed for max payload weight
vs2 = sqrt(W_NET/(0.5*rho_air*S*Cl_max3D)); % stall speed for no payload weight
vto2 = 1.2*vs2; % take off speed for no payload weight
