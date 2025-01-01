function diff_operators = getDiffOperators(k, sigma, N)
%	getDiffOperators       get differentiation operators
%
%       J. Mac Huang, 06/18/2016
%       final version 04/30/2023
%       AML DT, sculpt candy project
%
%   Helmholtz equation: u_xx + u_yy - sigma*u = f(x,y)
%
%   Inputs:
%            k            - order of IBSE method
%            sigma        - sigma in eqn above
%            N            - number of grids along one direction
%
%   Outputs: diff_operators containing the following data
%            Laplacian, invLaplacian - Laplacian and its inverse
%            Helmholtz_Op, inv_Helmholtz_Op
%                                    - Helmholtz Ops and its inverse
%            inv_H                   - inverse of H
%

% generate 2D Laplacian operator and its inverse in spectral space
j1 = repmat([0:N/2 -(N/2-1):-1],N,1);
j2 = repmat([0:N/2 -(N/2-1):-1]',1,N);
Laplacian = - (j1.^2 + j2.^2);
invLaplacian = 1./(Laplacian); invLaplacian(1,1) = 0;

dx = 1i*j1; dy = 1i*j2; 
inv_dx = 1./dx; inv_dx(:,1) = 0; 
inv_dy = 1./dy; inv_dy(1,:) = 0; 

% generate 2D Helmholtz operator and its inverse in spectral space
if sigma == 0
    Helmholtz_Op = Laplacian;
    inv_Helmholtz_Op = invLaplacian;
else
    Helmholtz_Op = (1-Laplacian*sigma);
    inv_Helmholtz_Op = 1./(Helmholtz_Op);
end

num_grid = 5; 
% generate the inverse of H operator for smooth extension
% Theta = max(1, 0.001*2^(-52)*(N/2)^(2*k+2));
%Theta = max(1, 1/(100*(2*pi/N))^(2*k+2));
Theta = (1/(num_grid*2*pi/N))^(2*(k+1));
inv_H = 1./(Laplacian.^(k+1) +Theta*(-1)^(k+1));

%Theta = max(1, 1/(100*(2*pi/N))^(2*k));
Theta = (1/(num_grid*2*pi/N))^(2*(k));
inv_H_lower = 1./(Laplacian.^(k) + Theta*(-1)^(k));

Theta = (1/(num_grid*2*pi/N))^(2*(k+2));
inv_H_higher = 1./(Laplacian.^(k+2) + Theta*(-1)^(k+2));

diff_operators = {Laplacian, invLaplacian, Helmholtz_Op, inv_Helmholtz_Op, inv_H, N, sigma, k, dx, dy, inv_dx, inv_dy, inv_H_lower, inv_H_higher};
end
