function ff = advective(N, u, v, f)
%	advective   Compute the advection term with exp(-36) anti-aliasing rule 
%
%       J. Mac Huang, 12/23/2015
%       final version 02/18/2016
%       AML DT, thermal transistor project
%
%   advective term: ff = u*d_x*f + v*d_y*f, see P33 of Peyret
%
%   Inputs:  N            - polynomial degree in x and y direction, must be
%                           odd
%            f            - function being advected
%            sigma        - constant
%
%   Outputs: ff           - advection term
%
%

f = reshape(f, N, N); u = reshape(u, N, N); v = reshape(v, N, N);

% get partial derivatives
ff = fft2(f); 
j1 = repmat([0:N/2 -(N/2-1):-1],N,1);
j2 = repmat([0:N/2 -(N/2-1):-1]',1,N);
    


% take variables into Fourier (Chebyshev) space
% ubar = fft2(u); vbar = fft2(v); 
fxbar = 1i*j1.*exp(-36*(j1./(N/2)).^36).*ff;
fybar = 1i*j2.*exp(-36*(j2./(N/2)).^36).*ff;

% gv = exp(-0.5*(jj1.^2+jj2.^2)*sigma^2);

ff = u.*(real(ifft2(fxbar))) + v.*(real(ifft2(fybar)));

ff = ff(:);
end