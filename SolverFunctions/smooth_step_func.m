function h = smooth_step_func(x,eps)
%	smooth_step_func       C^3 step function with Wendland functions
%
%       J. Mac Huang, 09/15/2022
%       AML DT, Plate Tectonics Project
%
%
%   Inputs:  x            - coordinate
%            eps          - half-width of the transition zone.

%
%   Outputs: h            - regularized Heaviside function
%

r = x/eps;
h = -3/4*r.^5 + 5/2*r.^4.*sign(r) - 5/2.*r.^3 + 5/4*r+1/2;
h = h.*(abs(r)<=1) + (r>1);

end