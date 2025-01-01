function [u, v, p] = VPSolver(boundary_data, boundary_condition_uv, diff_operators, operators, fu, fv, fp, invGrandSC)
%	VPSolver  Solving Stokes problem in the physical domain Omega
%
%       J. Mac Huang, 06/22/2023
%       AML DT, sculpt candy project
%
%   Stokes Problem:   L*u + p_x = XO*fu
%                     L*v + p_y = XO*fv
%                     u_x + v_y = XO*fp
%
%                     L is the Helmholtz operator, Helmholtz_Op
%
%   Inputs:  boundary_data, diff_operators, operators, inverse of SC 
%            generated by corresponding functions
%            fu, fv, fp   - driving terms defined in the equations above
%            boundary_condition_uv
%                         - boundary condition on u and v
%            invGrandSC   - inverse of grand SC
%
%   Outputs: u, v         - flow velocity in x, y direction
%            p            - pressure (zero mean)
%
[u_gamma, v_gamma, u_w1, v_w1, u_w2, v_w2] = boundary_condition_uv{:};
u_BC = [u_gamma; u_w1; u_w2]; v_BC = [v_gamma; v_w1; v_w2];

% grab all the data
[X, Y, nx, ny, Xw1, Yw1, Xw2, Yw2, XE, XEw, XO, ...
    theta, L, L, x0, y0, x, y, gamma] = boundary_data{:};
[Laplacian, invLaplacian, Helmholtz_Op, inv_Helmholtz_Op, inv_H, N, sigma,...
    IBSEk, dx, dy, inv_dx, inv_dy, inv_H_lower, inv_H_higher] ...
                    = diff_operators{:};
[Sn, normD, S, ST, STn, T, TT, T_1, T_2, T_3, R_1, R_2, R_3, B] ...
                    = operators{:};
                
% get all the length                
nbdy = length(X); nwall = length(Xw1); N = sqrt(length(XO));

% detect the order of method
if  IBSEk ==1
    Tn = T_1; R = R_1; Tn_lower = S; TT_lower = ST;
elseif IBSEk ==2
    Tn = T_2; R = R_2; Tn_lower = T_1; TT_lower = [ST; R_1];
elseif IBSEk == 3
    Tn = T_3; R = R_3; Tn_lower = T_2; TT_lower = [ST; R_2];
else 
    error('unknown order of IBSE')
end

% get other vector length              
length_T_u = size(Tn,2); 

% solve 0-Stokes problem (36a-c)              
[u0, v0, p0] = StokesSolver(boundary_data, diff_operators, ...
                            operators, XO.*fu, XO.*fv, XO.*fp);

% compute the RHS
RHS = [ST*u0-u_BC; ST*v0-v_BC; R*u0; R*v0; TT_lower*p0];

% get singular force
F = invGrandSC*RHS; 
Fu = F(1:length_T_u); Fv = F(length_T_u+1:2*length_T_u); 
Fp = F(2*length_T_u+1:end); 

% compute all xi values (25c-d)
xi_u = operator(-Tn*Fu, inv_H); xi_v = operator(-Tn*Fv, inv_H);
xi_p = operator(-Tn_lower*Fp, inv_H_lower); 

% finally, compute u, v, p from (25a-b)
fp_new = (1-XO).*(operator(xi_u, dx)+operator(xi_v, dy)) + XO.*fp;
fu_new = (1-XO).*(operator(xi_u, Helmholtz_Op)+operator(xi_p, dx)) + XO.*fu;
fv_new = (1-XO).*(operator(xi_v, Helmholtz_Op)+operator(xi_p, dy)) + XO.*fv;
[u, v, p] = StokesSolver(boundary_data, diff_operators, operators, fu_new, fv_new, fp_new);
end
