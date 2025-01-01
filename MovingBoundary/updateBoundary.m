function [boundary_data, V_data, renewflag] = ...
    updateBoundary(boundary_data, V_data, Vn, a_ice, epsilon, dt)
%   updateBoundary  update the boundary shape
%
%       J. Mac Huang, 12/23/2015
%       final version 02/18/2016
%       AML DT, thermal transistor project
%
%
%
%   Inputs: s       - dimensionless arclength
%           theta   - tangent angle
%           L       - total arclength
%           flux    - flux
%           beta    - dissolving speed coefficient (inverse latent heat)
%           epsilon - Gibbs Thompson effect coefficient
%
%   Output: Vtheta  - velocity for tangent angle
%           VL      - changing rate of arclength
%           Vy      - velocity of the top boundary point
%           implicitVtheta      - the diffusion part of theta equation
%           theta0  - periodic version of theta
%
%

% grab all the data
[X, Y, nx, ny, XE, XO,...
    theta, L, L_old, x0, y0, x, y, gamma, ~] = boundary_data{:};
[Vtheta, VL, Vn_anchor, a_ice_old, u_ice_old] = V_data{:};

% do time stepping
Vtheta_old = Vtheta; VL_old = VL; 
Vx_old = Vn_anchor*nx(1)+u_ice_old;  Vy_old = Vn_anchor*ny(1);
L_oldold= L_old; L_old = L;

% get all the vector length
nbdy = length(X);  N = sqrt(length(XO));

% define arclength coordinate and mean-free theta value
s= (0:nbdy-1)'/nbdy;
theta_mf = theta+2*pi*s;

% determine the grid point number for the new boundary
nbdyold = nbdy;
nbdy = floor(L/(2*pi/N)/2); if mod(nbdy,2)~=0; nbdy = nbdy+1;end

% interpolate old data to new boundary
s= (0:nbdy-1)'/nbdy;
theta_mf  = interpft(theta_mf, nbdy);
Vtheta_old = interpft(Vtheta_old, nbdy);
Vn = interpft(Vn, nbdy);

% if boundary size has changed, then renew the SC
if nbdy~=nbdyold
    renewflag = 1;
    disp( '  boundary size changed, SC will be renewed')
else
    renewflag = 0;
end

% computing all the boundary velocity through theta-L method
j1 = [0:nbdy/2 -(nbdy/2-1):-1]'*2*pi;
dtheta_ds = real(ifft(1i*j1.*fft(theta_mf)))-2*pi;
Vn_tmp = Vn + epsilon/L*(dtheta_ds + 2*pi);
[intVn, Vnend] = fourierInt((dtheta_ds).*Vn_tmp);
Vs = intVn - s*Vnend;

Vtheta = 1/L*(real(ifft(1i*j1.*fft(Vn)))+(dtheta_ds).*Vs);
VL = -Vnend;
Vn_anchor = -Vn_tmp(1);
u_ice = u_ice_old +  dt/2*(3*a_ice - a_ice_old);
Vx = Vn_anchor*nx(1)+u_ice;  Vy = Vn_anchor*ny(1);

% evolve the boundary
L = L + dt/2*(3*VL - VL_old);
x0 = x0 + dt/2*(3*Vx-Vx_old);
y0 = y0 + dt/2*(3*Vy-Vy_old);
gSigma1 = 2*pi*sqrt(epsilon*(dt))*sqrt(L.^(-3/2)+L_old.^(-3/2));
gSigma2 = 2*pi*sqrt(epsilon*(dt))*sqrt(L.^(-3/2)+2*L_old.^(-3/2) +L_oldold.^(-3/2));
N1 =gaussianF( Vtheta, gSigma1) ; N2 = gaussianF(Vtheta_old, gSigma2);
theta_mf = gaussianF( theta_mf, gSigma1 ) + dt/2*(3*N1-N2);
theta = theta_mf - 2*pi*s;

% gegenerate the boundary data, keep a record of L_old
boundary_data = boundaryinfo(theta, L, x0, y0, x, y, gamma);
boundary_data{9} = L_old;

% save all the boundary velocity data
V_data = {Vtheta, VL, Vn_anchor, a_ice, u_ice, Vn};
end
