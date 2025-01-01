%parameter_initialize

% exp_parameter
inch = 2.54;
height_ratio = 1;                % normalize the height to be the ratio
H = 0.01*3*inch*sin(alpha)/height_ratio;      % ice height
W = 0.06;                          % ice width
T0 = 0;                            % melting temp
T_star = 4;                        % highest density temp
rho_star = 999.84;                 % density at 4
kt = 1.3e-7;                       % thermal diffusivity
tau = H^2/kt;                      % diffusivity time scale
k_viscosity = 1.57e-6;             % nu, kinematic viscosity
gravity = 9.8;                     % gravitational acceleration
beta = 7.68e-6;                    % curvature for rho(T) 
Cp = 4.205e3;                      % heat capacity
Latent = 3.33e5;                   % latent heat
rho_ice = 917;                     % density of ice at 0

% Dimensionless numbers
Pr = k_viscosity/kt;                                     % Prandtl 
Ra = gravity*beta*(T_infty-T0)^2*H^3/(k_viscosity*kt);   % Rayleigh 
St = Cp*(T_infty-T0)/Latent;                             % Stefan
lambda = rho_ice/rho_star;                           % dimensionless density
lambda = 9;
lambda = rho_ice/rho_star*10;

% Solver parameters
gamma = 4/3;                      % aspect ratio of physical domain
q = 0.0;                          % concentration depletion coeff.
endFraction = 0.1;                % smallest possible size (for shrinking)
gmres_tol = 1e-6;                 % tolerance for the GMRES
iter_max = 10;                    % max GMRES iterations for SC to renew
epsilon = .03;                   % Gibbs Thomson effect
damping_strength = 0;             % how strong the damping is
damping_size = 0.8;               % how wide the damping is

% damping mask
width = 1;
y_lower = width*1.2; y_upper = 2*pi - width*1.2;
damp_strength = 400;

