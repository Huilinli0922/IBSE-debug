function [u, iter_u] = GMRES_HelmholtzSolver(boundary_data, ...
            boundary_condition, diff_operators, operators, f, SC, gmres_tol)
%	GMRES_HelmholtzSolver  Helmholtz solver in a channel,
%                           with GMRES and previous SC preconditioner.
%
%       J. Mac Huang, 06/18/2016
%       final version 06/19/2016
%       AML DT, sculpt candy project
%
%   Helmholtz equation:  u - sigma*(u_xx + u_yy) = f(x,y)
%
%   Inputs:  boundary_data, diff_operators, operators, inverse of SC 
%            generated by corresponding functions
%            f            - RHS of the Helmholtz eqn
%            SC           - Schur complement
%            gmres_tol    - tolerance for gmres
%
%
%   Outputs: u            - solution at (X,Y)
%            iter_u       - iterations for GMRES to complete
%
%

% grab all the data
[X, Y, nx, ny, Xw1, Yw1, Xw2, Yw2, XE, XEw, XO] = boundary_data{:};
[Laplacian, invLaplacian, Helmholtz_Op, ...
    inv_Helmholtz_Op, inv_H, N, sigma, IBSEk] = diff_operators{:};
[Sn, normD, S, ST, STn, T, TT, T_1, T_2, T_3, R_1, R_2, R_3, B] ...
    = operators{:};
[a0, b0, g0, a_w1, b_w1, g_w1, a_w2, b_w2, g_w2] = boundary_condition{:};
g = [g0; g_w1; g_w2];


if IBSEk == 1
    Tn = T_1; R = R_1;
elseif IBSEk == 2
    Tn = T_2; R = R_2;
elseif IBSEk == 3
    Tn = T_3; R = R_3;
else
    error('unknown order of IBSE')
end

% compute the RHS, so SC*F = RHS
XE = 1-XO;
invLf = operator( (XO).*f, inv_Helmholtz_Op);
RHS = [ -R*(invLf); B*(invLf) - g];

% compute singular force F
[F, ~, ~, iter_u] = gmres(@forward_op, RHS, 20, gmres_tol, 100, SC);

% from singular force F, compute xi and u
xi = -operator(Tn*F, inv_H);
u = (operator(f.*(XO) + (1-XO).*operator(xi, Helmholtz_Op),inv_Helmholtz_Op));

% the function below is to form the operation SC forwardly, for using GMRES
    function newRHST = forward_op(F)
        xi_tmp = operator(Tn*F, inv_H);
        u_tmp = operator( XE.*operator(xi_tmp, Helmholtz_Op), inv_Helmholtz_Op);
        newRHST = [R*(xi_tmp-u_tmp); B*u_tmp];
    end
end



