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

[X, Y, nx, ny, Xw1, Yw1, Xw2, Yw2, XE, XEw, XO] = boundary_data{:};
% gamma length and wall length
nbdy = length(X); nwall = length(Xw1);

% preallot space for all the operators
S = sparse(N^2, length(X)); 
T0 = S; T1 = S; T2 = S; T3 = S; TT0 = T0; TT1 = T1; TT2 = T2; TT3 = T3;
S_w1 = sparse(N^2, length(Xw1)); 
T0_w1 = S_w1; T1_w1 = S_w1; T2_w1 = S_w1; T3_w1 = S_w1;
TT0_w1 = S_w1; TT1_w1 = S_w1; TT2_w1 = S_w1; TT3_w1 = S_w1;
S_w2 = sparse(N^2, length(Xw2)); 
T0_w2 = S_w2; T1_w2 = S_w2; T2_w2 = S_w2; T3_w2 = S_w2;
TT0_w2 = S_w2; TT1_w2 = S_w2; TT2_w2 = S_w2; TT3_w2 = S_w2; 

% get all the regularized delta functions
delta_functions = delta3(N, X, Y);
delta_functions_wall_1 = delta3(N, Xw1, Yw1);
delta_functions_wall_2 = delta3(N, Xw2, Yw2);

[d0x, d1x, d2x, d3x, d0y, d1y, d2y, d3y] = delta_functions{:};
[d0xw1, d1xw1, d2xw1, d3xw1, d0yw1, d1yw1, d2yw1, d3yw1] = delta_functions_wall_1{:};
[d0xw2, d1xw2, d2xw2, d3xw2, d0yw2, d1yw2, d2yw2, d3yw2] = delta_functions_wall_2{:};


% grab boundary conditions
[a0, b0, g0, a_w1, b_w1, g_w1, a_w2, b_w2, g_w2] = boundary_condition{:};


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

% on bottom wall
for k = 1:length(Xw1)
    
    % get all the products of delta functions
    d0xd0y_w1 = kron(d0xw1(k,:), d0yw1(k,:)');
    d0xd1y_w1 = kron(d0xw1(k,:), d1yw1(k,:)');
    d0xd2y_w1 = kron(d0xw1(k,:), d2yw1(k,:)');
    d0xd3y_w1 = kron(d0xw1(k,:), d3yw1(k,:)');
    
    % get all the operators
    S_w1(:,k) = d0xd0y_w1(:);
    T0_w1(:,k) = S_w1(:,k);
    T1_w1(:,k) = -d0xd1y_w1(:);
    T2_w1(:,k) = d0xd2y_w1(:);  
    T3_w1(:,k) = -d0xd3y_w1(:);  
    TT0_w1(:,k) = T0_w1(:,k);
    TT1_w1(:,k) = T1_w1(:,k);
    TT2_w1(:,k) = T2_w1(:,k);
    TT3_w1(:,k) = T3_w1(:,k);
end

% on top wall
for k = 1:length(Xw2)
    
    % get all the products of delta functions
    d0xd0y_w2 = kron(d0xw2(k,:), d0yw2(k,:)');
    d0xd1y_w2 = kron(d0xw2(k,:), d1yw2(k,:)');
    d0xd2y_w2 = kron(d0xw2(k,:), d2yw2(k,:)');
    d0xd3y_w2 = kron(d0xw2(k,:), d3yw2(k,:)');
    
    % get all the operators
    S_w2(:,k) = d0xd0y_w2(:);
    T0_w2(:,k) = S_w2(:,k);
    T1_w2(:,k) = -d0xd1y_w2(:);
    T2_w2(:,k) = d0xd2y_w2(:);
    T3_w2(:,k) = -d0xd3y_w2(:);
    TT0_w2(:,k) = T0_w2(:,k);
    TT1_w2(:,k) = T1_w2(:,k);
    TT2_w2(:,k) = T2_w2(:,k);
    TT3_w2(:,k) = T3_w2(:,k);
end


% put all operators together
T_w1 = [T0_w1 T1_w1 T2_w1 T3_w1];                % T on bottom wall
TT_w1 = [TT0_w1 TT1_w1 TT2_w1 TT3_w1]';          % T^* on bottom wall
ST_w1 = TT0_w1';                              % S^* on bottom wall

T_w2 = [T0_w2 T1_w2 T2_w2 T3_w2];                % T on top wall
TT_w2 = [TT0_w2 TT1_w2 TT2_w2 TT3_w2]';          % T^* on top wall
ST_w2 = TT0_w2';                              % S^* on top wall
    
T = [T0 T1 T2 T3];                          % T on Gamma
ST = TT0';                                  % S^* on Gamma
TT = [TT0 TT1 TT2 TT3]';                    % T^* on Gamma

S = [S S_w1 S_w2]; ST = [ST; ST_w1; ST_w2];     % put S, S^* together
T = [T T_w1 T_w2]; TT = [TT; TT_w1; TT_w2];     % put T, T^* together


% Sn and STn are the S and S^* for Dirichlet on Gamma and Neumann on walls,
%             this is the case for the thermal and mass diffusion study
Sn = sparse([S T1_w1 T1_w2]); STn = sparse([ST; TT1_w1'; TT1_w2']);


% normal derivative operator, something useful
normD = sparse([TT1'; TT1_w1'; TT1_w2']);

% get boundary condtion operators
B_0 = a0*TT0'+ b0*TT1'; 
B_w1 = a_w1*TT0_w1' + b_w1*TT1_w1'; 
B_w2 = a_w2*TT0_w2' + b_w2*TT1_w2';
B = sparse([B_0; B_w1; B_w2]);


% T_k and R_k are T and R for IBSE k method
T_1 = T(:,[1:2*nbdy 4*nbdy+1:4*nbdy+2*nwall  4*nbdy+4*nwall+1:4*nbdy+6*nwall]);
T_2 = T(:,[1:3*nbdy 4*nbdy+1:4*nbdy+3*nwall  4*nbdy+4*nwall+1:4*nbdy+7*nwall]);
T_3 = T;
R_1 = TT([nbdy+1:nbdy*2 (4*nbdy+nwall+1):(4*nbdy+2*nwall) (4*nbdy+5*nwall+1):(4*nbdy+6*nwall)], :);
R_2 = TT([nbdy+1:nbdy*3 (4*nbdy+nwall+1):(4*nbdy+3*nwall) (4*nbdy+5*nwall+1):(4*nbdy+7*nwall)], :);
R_3 = TT([nbdy+1:nbdy*4 (4*nbdy+nwall+1):(4*nbdy+4*nwall) (4*nbdy+5*nwall+1):(4*nbdy+8*nwall)], :);

% collect all the operators
operators = {Sn, normD, S, ST, STn, T, TT, T_1, T_2, T_3, R_1, R_2, R_3, B};
end
