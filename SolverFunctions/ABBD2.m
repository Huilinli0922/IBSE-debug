function F = ABBD2(Dt,f, B, Bold, omega_old,omega_oldold, ifstart)
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
%           ifstart      - run Euler if this is the first step
%
%   Output: F - RHS for the Helmholtz equation
%
%

if ifstart
    F = Dt*(f - B) + omega_old;
else
    F = 2*Dt/3*(f - (2*B-Bold) + (4*omega_old-omega_oldold)/(2*Dt));
end

end
