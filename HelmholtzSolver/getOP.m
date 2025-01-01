function operators = getOP(boundary_data, boundary_condition, N)
%	getOP       get operators for the Helmholtz solver
%
%       J. Mac Huang, 06/18/2016
%       final version 04/30/2023
%       AML DT, sculpt candy project
%
%   Helmholtz equation: u_xx + u_yy - sigma*u = f(x,y)
%
%   Inputs:  
%            N            - grid size
%
%            boundary_data contains the following variables            
%            X, Y         - Gamma boundary points
%            nx, ny       - normal direction of the boundary
%            ds           - metric on the boundary, ds = || dX/d(theta) ||
%            Xw1, Yw1     - nodes on the lower wall
%            Xw2, Yw2     - nodes on the upper wall
%
%   Outputs: operators containing the following data
%            S, ST, T, TT - spread and interpolation operators
%                           for T and TT the order of element is 
%                           [Gamma boundary, lower wall, upper wall]
%            STn          - interpolation operator for Dirichlet cond. on
%                           moving boundary and Neumann cond. on wall
%            normD        - take gradient at the boundary
%

[X, Y, nx, ny, XE, XO] = boundary_data{:};
% gamma length and wall length
nbdy = length(X); 

% preallot space for all the operators
S = sparse(N^2, length(X)); 
T0 = S; T1 = S; T2 = S; T3 = S; TT0 = T0; TT1 = T1; TT2 = T2; TT3 = T3;


% get all the regularized delta functions
delta_functions = delta3(N, X, Y);


[d0x, d1x, d2x, d3x, d0y, d1y, d2y, d3y] = delta_functions{:};


% grab boundary conditions
[a0, b0, g0] = boundary_condition{:};


nxy = nx.*ny; nxx = nx.*nx; nyy = ny.*ny;
nxxx = nxx.*nx; nxxy = nxx.*ny; nyyy = nyy.*ny; nxyy = nxy.*ny;


% on boundary Gamma
for k = 1:length(X)
    
    % get all the products of delta functions
    d0xd0y = kron(d0x(k,:), d0y(k,:)');
    d1xd0y = kron(d1x(k,:), d0y(k,:)');
    d0xd1y = kron(d0x(k,:), d1y(k,:)');
    d2xd0y = kron(d2x(k,:), d0y(k,:)');
    d0xd2y = kron(d0x(k,:), d2y(k,:)');
    d1xd1y = kron(d1x(k,:), d1y(k,:)');
    d3xd0y = kron(d3x(k,:), d0y(k,:)');
    d2xd1y = kron(d2x(k,:), d1y(k,:)');
    d1xd2y = kron(d1x(k,:), d2y(k,:)');
    d0xd3y = kron(d0x(k,:), d3y(k,:)');

    % get all the operators
    S(:,k) = d0xd0y(:);
    T0(:,k) = S(:,k);
    T1(:,k) = - d1xd0y(:)*nx(k) - d0xd1y(:)*ny(k);
    T2(:,k) =   d2xd0y(:)*nxx(k) + d0xd2y(:)*nyy(k) + 2*d1xd1y(:)*nxy(k);
    T3(:,k) = - d3xd0y(:)*nxxx(k) - d0xd3y(:)*nyyy(k) - 3*d1xd2y(:)*nxyy(k) ...
                 - 3*d2xd1y(:)*nxxy(k);
    TT0(:,k) = S(:,k);
    TT1(:,k) = T1(:,k);
    TT2(:,k) = T2(:,k);
    TT3(:,k) = T3(:,k);
end

% put all operators together

T = [T0 T1 T2 T3];                          % T on Gamma
ST = TT0';                                  % S^* on Gamma
TT = [TT0 TT1 TT2 TT3]';                    % T^* on Gamma


% Sn and STn are the S and S^* for Dirichlet on Gamma and Neumann on walls,
%             this is the case for the thermal and mass diffusion study
Sn = sparse(S); STn = sparse(ST);


% normal derivative operator, something useful
normD = sparse(TT1');

% get boundary condtion operators
B_0 = a0*TT0'+ b0*TT1'; 
B = sparse(B_0);


% T_k and R_k are T and R for IBSE k method
T_1 = T(:,1:2*nbdy);
T_2 = T(:,1:3*nbdy);
T_3 = T;
R_1 = TT(nbdy+1:nbdy*2, :);
R_2 = TT(nbdy+1:nbdy*3, :);
R_3 = TT(nbdy+1:nbdy*4, :);

% collect all the operators
operators = {Sn, normD, S, ST, STn, T, TT, T_1, T_2, T_3, R_1, R_2, R_3, B};
end
