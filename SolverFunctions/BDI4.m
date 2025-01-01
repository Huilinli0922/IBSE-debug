function F = BDI4(omega_old, omega_oldold, omega_oldoldold, omega_oldoldoldold)
%   ABBD2   Second order Adams-Bashforth/backward-difference method
%   (this is a function modified from Zoltn Csti's code)
%
%       J. Mac Huang, 12/23/2015
%       final version 02/18/2016
%       AML DT, thermal transistor project
%
%   advective term: ff = u*d_x*f + v*d_y*f
%
%   Inputs: Dt    - time step
%           f     - forcing term
%           B     - convective term
%           Bold  - convective term one time-step before
%           omega_old    - vorticity one time-step before
%           omega_oldold - vorticity two time-steps before
%           nu    - kinematic viscosity
%
%   Output: F - RHS for the Helmholtz equation
%
%

F = 12/25*(4*omega_old-3*omega_oldold+4/3*omega_oldoldold-1/4*omega_oldoldoldold);
end
