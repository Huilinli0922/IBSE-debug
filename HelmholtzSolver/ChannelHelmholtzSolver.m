function u = ChannelHelmholtzSolver(boundary_data, boundary_condition, diff_operators, operators, f, invSC)
%	ChannelHelmholtzSolver  Helmholtz solver in a channel,
%                           need to run ChannelHelmholtzSolverPre first.
%
%       J. Mac Huang, 06/18/2016
%       final version 06/19/2016
%       AML DT, sculpt candy project
%
%   Helmholtz equation: u - sigma*(u_xx + u_yy) = f(x,y)
%
%   Inputs:  boundary_data, diff_operators, operators, inverse of SC 
%            generated by corresponding functions
%            F            - RHS of the Helmholtz eqn
%            g            - boundary condition
%            invSC        - inverse of SC
%
%
%   Outputs: u            - solution at (X,Y)
%


% grab all the data
[X, Y, nx, ny, Xw1, Yw1, Xw2, Yw2, XE, XEw, XO] = boundary_data{:};
[Laplacian, invLaplacian, Helmholtz_Op, inv_Helmholtz_Op, inv_H, N, sigma, IBSEk] ...
                    = diff_operators{:};
[Sn, normD, S, ST, STn, T, TT, T_1, T_2, T_3, R_1, R_2, R_3, B] ...
                    = operators{:};
[a0, b0, g0, a_w1, b_w1, g_w1, a_w2, b_w2, g_w2] = boundary_condition{:};
g = [g0; g_w1; g_w2];

nbdy = length(X); nwall = length(Xw1);     
if IBSEk == 1
    Tn = T_1; R = R_1;
elseif IBSEk == 2
    Tn = T_2; R = R_2;
elseif IBSEk == 3
    Tn = T_3; R = R_3;
else 
    error('unknown order of IBSE')
end

% from SC, compute F, and the solution
invLf = operator( (XO).*f, inv_Helmholtz_Op); 
RHS = [ R*(invLf); B*(invLf) - g];

if sigma == 0
    FG = invSC*[RHS; sum((XO).*f)/(nbdy+2*nwall)^2];
    F = FG(1:end-1); U = FG(end);
    xi = -operator(Tn*F, inv_H);
    u = (operator(f.*(XO) + (1-XO).*operator(xi, Helmholtz_Op),inv_Helmholtz_Op))-U/(nbdy+2*nwall);
else
    F = invSC*RHS;
    xi = -operator(Tn*F, inv_H);
    u = (operator(f.*(XO) + (1-XO).*operator(xi, Helmholtz_Op),inv_Helmholtz_Op));
end

end
