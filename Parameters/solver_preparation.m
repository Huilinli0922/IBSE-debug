% solver_preparation

% Initial shape
% initialize boundary geometry
x0 = pi; y0 = pi+1/2;                     % physical location of s=0
l0 = 1/sin(alpha); L0 = (l0+1+1/tan(alpha))*height_ratio;

% compute the number of grid points for candy boat
nbdy = floor(L0/(2*pi/N)/2);
if mod(nbdy,2)~=0; nbdy = nbdy+1; end

% define arclength and geometry
s = (0:nbdy-1)'/nbdy;


[theta, ~, ~, L] = closed_wedge(s, alpha, l0, 0.04,height_ratio);


% define boundary conditions
g = zeros(nbdy,1); 
c_wall=(T_infty-T_star)/(T_infty-T0);
c_ice=(T0-T_star)/(T_infty-T0);

boundary_condition_T = {1, 0, g+c_ice};
boundary_condition_uv = {g, g};

% define cartesian coordinates for the domain and boundaries
[x,y] = meshgrid((1:N)/N*2*pi);
x = x(:); y = y(:);


% define damping mask
[f, y_lower, y_upper] = damping_mask(y, y_lower, y_upper, width);
f = damp_strength*f;

%% Preprocessing stage
% define sigma in the Helmholtz eqns
nu = Pr; %1/Re;
diffusivity = 1; %1/Pe;
sigma_u = (2/3*dt*nu);
sigma_T = (2/3*dt*diffusivity);

% generate all boundary data
boundary_data = boundaryinfo(theta, L, x0, y0, x, y, gamma);
[X, Y, nx, ny, XE, XO, S0] = boundary_data{:};
m_ice = lambda*S0;

% differential operators
diff_operators_uv = getDiffOperators(1, sigma_u, N);
diff_operators_T = getDiffOperators(1, sigma_T, N);
[Laplacian, invLaplacian, Helmholtz_Op, inv_Helmholtz_Op, inv_H, N, sigma, ...
    k, dx, dy, inv_dx, inv_dy, inv_H_lower] = diff_operators_uv{:};

% initialize all variables
u = 0*x; v = 0*x; p = 0*x; c = 0*x+c_wall;
B_u = advective(N, u, v, u); B_v = advective(N, u, v, v);
C = advective(N, u, v, c);

cold = c; coldold = cold;
uold = u; uoldold = uold;
vold = v; voldold = vold;
Cold = C; Coldold = Cold;
B_u_old = B_u; B_u_oldold = B_u_old;
B_v_old = B_v; B_v_oldold = B_v_old;


% initialize all shape related variables
nbdy0 = nbdy; L0 = L;
V_data = {0*theta, 0, 0, 0, 0}; Vn = 0*theta; a_ice = 0; u_ice = 0;
% u_bd= 0; v_bd= 0; Vx= 0; Vy= 0;
renewflag1 = 1; renewflag2 = 1;

% prepare solver behaviors
pp_plot = 1; k = 1; pp_save = 1; save_num = 1; NANflag = 0;
pp_save2 = 1; save_num2 =1;
tt = 0; ifstart = 0;

% prepare all operators
operators = getOP(boundary_data, boundary_condition_T, N);
[Sn, normD, S, ST, STn, T, TT, T_1, T_2, T_3, R_1, R_2, R_3, B] ...
    = operators{:};

% prepare Shur Complement
disp('    preparing Shur complement...')