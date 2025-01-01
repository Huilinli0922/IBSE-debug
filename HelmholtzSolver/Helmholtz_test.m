
% this is a test code for the Helmholtz solver
clear all
close all

% test Helmholtz: u - sigma*(u_xx + u_yy) = f(x,y)

% parameters
N = 128;
nbdy = N/2;
nwall = N/2;
gamma = 1.2;
sigma = 1e-3;
k =1;


% generate boundary Gamma and meshgrid
s = (0:nbdy-1)'/nbdy; theta = pi/2-2*pi*s; L = 2*pi; x0 = pi; y0 = pi+1;
[x,y] = meshgrid((1:N)/N*2*pi);
x = x(:); y = y(:);
boundary_data = boundaryinfo(theta, L, x0, y0, x, y, gamma, nwall);
[X, Y, nx, ny, Xw1, Yw1, Xw2, Yw2, XE, XEw, XO] = boundary_data{:};


% exact solution we chose (make it periodic, nicer this way)
u_exact = 6*cos(2*y).*sin(x);


% associated RHS and boundary conditions
f = (1+5*sigma)*u_exact;

Dirichlet_Gamma = 6*cos(2*Y).*sin(X);
Dirichlet_w1 = 6*cos(2*Yw1).*sin(Xw1);
Dirichlet_w2 = 6*cos(2*Yw2).*sin(Xw2);

Neumann_Gamma = 6*cos(2*Y).*cos(X).*nx - 12*sin(2*Y).*sin(X).*ny;
Neumann_w1 = -12*sin(2*Yw1).*sin(Xw1);
Neumann_w2 = -12*sin(2*Yw2).*sin(Xw2);

a0 = 1; b0 = 0; 
g0 = a0*Dirichlet_Gamma + b0*Neumann_Gamma;

a_w1 = 1; b_w1 = 0; 
g_w1 = a_w1*Dirichlet_w1 + b_w1*Neumann_w1;

a_w2 = 1; b_w2 = 0; 
g_w2 = a_w2*Dirichlet_w2 + b_w2*Neumann_w2;

boundary_condition = {a0, b0, g0, a_w1, b_w1, g_w1, a_w2, b_w2, g_w2};


% get all the operators
diff_operators = getDiffOperators(k, sigma, N);
operators = getOP(boundary_data, boundary_condition, N);


% preprocessing, get SC
SCLU = ChannelHelmholtzSolverPre(boundary_data, diff_operators, operators);
inv_SCLU=inv(SCLU);

% compute
u = ChannelHelmholtzSolver(boundary_data, boundary_condition, diff_operators, operators, f, inv_SCLU);


% show the solution
hh = pcolor(reshape(x,N,N), reshape(y,N,N), reshape(u.*XO,N,N));

set(hh, 'EdgeColor','none')
shading interp
colorbar; 
axis equal
axis off

hold on
plot(Xw1, Yw1, 'r', 'lineWidth', 2)
plot(Xw2, Yw2, 'r', 'lineWidth', 2)
plot([X; X(1)], [Y; Y(1)], 'r', 'lineWidth', 2)
hold off

title('solution, u')


% show the error
figure
hh = pcolor(reshape(x,N,N), reshape(y,N,N), reshape(abs(u-u_exact).*XO,N,N));
set(hh, 'EdgeColor','none')
shading interp
colorbar;
axis equal
axis off

hold on
plot(Xw1, Yw1, 'r', 'lineWidth', 2)
plot(Xw2, Yw2, 'r', 'lineWidth', 2)
plot([X; X(1)], [Y; Y(1)], 'r', 'lineWidth', 2)
hold off

title(['error, L^\infty error = ' num2str(norm((u-u_exact).*XO,'inf'))])