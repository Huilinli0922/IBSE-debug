%   ice boat VP simulate the ice melting problem with
%   velocity-pressure format
%   J. Mac Huang, 06/23/2023
%   AML DT, candy boat project
%   Note: this function requires parallel computing toolbox

clear all; close all; clear global all 

fftw('wisdom'); 

% Server behavior and initialize

location = 2;                           % 0 cims server, 1 E605, 2 Tiffany, 
                                        % 3 Mac home, 4 nyush_hpc, 5
                                        % nyush_tiff, 

N = 256;                                % node number in x and y direction
alpha = pi/6;                           % inclination angle
disp(['alpha=' num2str(round(alpha*180/pi))])
T_infty = 10;                           % ambient temp

save_skip = 200;                       % save every # time steps
disk_save_int = 1000*save_skip;%1000      % how many steps until disk save
disk_save_int2 = 1000; %100000
renew_step=[];
hold_step = 100;

% Add path, prepare everything
addpath('./SolverFunctions', './HelmholtzSolver', ...
    './VelocityPressureSolver', './VisualizationTools', ...
    './IceShape', './MovingBoundary', './Parameters');

server_plot_save_setup                  % prepare plot and save
parameter_initialize                    % other variable + exp parameter

dt = 0.001/sqrt(Ra);                              % time step
% dt = 0.002/sqrt(Ra); 

solver_preparation                      % prepare shape and operators



%% main loop
while L/L0 > endFraction

    [boundary_data, V_data, renewflag1] = ...
        updateBoundary(boundary_data, V_data, Vn, a_ice, epsilon, dt);
    [X, Y, nx, ny, XE, XO,...
        theta, L, ~, x0, y0, ~, ~, ~, X_np,S0] = boundary_data{:};

    u_ice = V_data{5};                  % ice moving velocity
   
    if k < hold_step
	   Vn = V_data{6}*0;
   else
	   Vn = V_data{6};
   end

   % Vn = V_data{6};                     % boundary normal velocity
    m_ice = lambda*S0;                  % mass of ice

    % define the new boundary conditions
    nbdy = length(X); g = zeros(nbdy,1);
    u_BC = u_ice + Vn.*nx; 
    v_BC = Vn.*ny;
    boundary_condition_T = {1, 0, g+c_ice};
    boundary_condition_uv = {u_BC, v_BC};

    % update the oprators
    operators = getOP(boundary_data, boundary_condition_T, N);
    [Sn, normD, S, ST, STn, T, TT, T_1, T_2, T_3, R_1, R_2, R_3, B] = operators{:};

    % renew the SC if needed
    % flag1 - changed nbdy as boundary evloves; flag2 - too many GMRES iterations
    if renewflag1 || renewflag2
        renew_step =[renew_step, tt];
        disp('renewing SC...')
        % update GrandSC
        GrandSC = VPPre(boundary_data, diff_operators_uv, operators);
        % update SC
        SC = ChannelHelmholtzSolverPre(boundary_data, diff_operators_T, operators);
        % reset flags
        renewflag1 = 0; renewflag2 = 0;
    end

    % keep a record of variables in previous steps
    coldold = cold; cold = c; uoldold = uold; uold = u;
    voldold = vold; vold = v; B_u_old = B_u; B_v_old = B_v; Cold = C;

    % currently we are at time step n
    % c, u, v, cold, uold, vold  are still at (n-1) step
    % coldold, uoldold, voldold, B_*_old, Cold at (n-2) step

    % solve for the concentration/temperature field
    C = advective(N, u, v, c) + f.*(c - c_wall);   % C at (n-1) step
    RHST = ABBD2(dt,0, C, Cold, cold, coldold, ifstart);
    [c, iter_c] = GMRES_HelmholtzSolver(boundary_data, ...
        boundary_condition_T, diff_operators_T, operators, RHST, SC, gmres_tol);
    Vn = -St*normD*c; Vn = Vn(1:nbdy);
    % this is c at step n

    % solve for the flow speed (u, v) and pressure p
    B_u = advective(N, u, v, u) + f.*u;
    B_v = advective(N, u, v, v) + f.*v;       % now B_* are at step (n-1)

    fu = ABBD2(dt, 0, B_u, B_u_old, uold, uoldold, ifstart);
    fv = ABBD2(dt, Pr*Ra*c.^2, B_v, B_v_old, vold, voldold, ifstart);
    fp = 0;

    [u, v, p, iter_u] = GMRES_VPSolver(boundary_data, boundary_condition_uv, ...
        diff_operators_uv, operators, fu, fv, fp, GrandSC, gmres_tol);

    % get the pressure distribution on ice surface
    % note: pressure p is rescaled by 2*dt/3, due to the way ABBD2 works
    % p is also shifted so the top wall has 0 mean pressure

    % map p to boundary surf_p
    p = p/(2/3*dt); surf_p = ST*p;

    % % mean pressure on top wall
    % wall_p_up = mean(surf_p(end-nwall+1:end));
    % 
    % % mean pressure on buttom wall
    % wall_p_down = mean(surf_p(end-2*nwall+1:end-nwall));
    % 
    % % pressure gradient vertically
    % p_grad = (wall_p_down - wall_p_up)/(Yw2(1)-Yw1(1)); p_grad
    % 
    % % make surf_p gradient free and mean free
    % surf_p = surf_p(1:nbdy) + p_grad*Y; surf_p = surf_p - mean(surf_p);

    % calc stress tensor
    surf_tau = Pr*(ST*(operator(u,dy)+operator(v,dx)));
    surf_tau = surf_tau(1:nbdy);
    diag_stuff = 2*Pr*(ST*(operator(u,dx)));
    diag_stuff = diag_stuff(1:nbdy);

    pressure_field=[];
    shear_field=[];
    pressure_field=[-surf_p.*nx,-surf_p.*ny];
    shear_field=[diag_stuff.*nx+surf_tau.*ny,-diag_stuff.*ny+surf_tau.*nx];
    pressure=[mean(-surf_p.*nx), mean(-surf_p.*ny)];
    shear=[mean(diag_stuff.*nx+surf_tau.*ny), mean(-diag_stuff.*ny+surf_tau.*nx)];
    f_x=(pressure(1)+shear(1))*L;
    f_y=(pressure(2)+shear(2))*L;

    % get ice accelaration
    if k < hold_step
        a_ice = 0;
    else
        a_ice = f_x/m_ice;
    end

    % check nan or inf
    if isnan(norm(c,Inf))
        error('Did not converge')
        break
    end

    plot_print_save
    k = k+1; tt = tt+dt;

end
