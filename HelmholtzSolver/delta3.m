function delta_functions = delta3(N, X, Y)
%	delta3  generates discrete Delta function and its derivatives.
%           (x,y) in C, (X,Y) on boundary Gamma
%
%       J. Mac Huang, 06/18/2016
%       final version 04/30/2023
%       AML DT, sculpt candy project
%
%   Inputs:  x, y         - meshgrid solution at (X,Y)
%            X, Y         - boundary points
%            dx           - grid size, y(2) - y(1)
%
%   Outputs:
%            delta_functions contain the following variables
%            dix, diy     - i th derivative of the delta function for x, y
%


% make computational grids
x = (1:N)'/N*2*pi; dx = x(2) - x(1);

% as the domain is periodic, extend it by 10 grid points so delta functions
% can be computed on this extended domain. We will send everything back to
% the smaller domain later
xextend = [x(end-8:end) - 2*pi; x; x(1:9) + 2*pi];
yextend = xextend;

% the meshgrid of (x-X) and (y-Y)
rx = (repmat(xextend',length(X),1)-repmat(X,1,length(xextend)));
ry = (repmat(yextend',length(X),1)-repmat(Y,1,length(yextend)));

% get all the derivatives
d0x = stencil_4th_C3(rx,dx,0);
d1x = stencil_4th_C3(rx,dx,1);
d2x = stencil_4th_C3(rx,dx,2);
d3x = stencil_4th_C3(rx,dx,3);
d0y = stencil_4th_C3(ry,dx,0);
d1y = stencil_4th_C3(ry,dx,1);
d2y = stencil_4th_C3(ry,dx,2);
d3y = stencil_4th_C3(ry,dx,3);

% contert the extended domain back to the smaller one
d0x = d0x(:, 10:end-9) + [zeros(length(X), N-9) d0x(:, 1:9)] + [d0x(:, end-8:end) zeros(length(X), N-9) ];
d1x = d1x(:, 10:end-9) + [zeros(length(X), N-9) d1x(:, 1:9)] + [d1x(:, end-8:end) zeros(length(X), N-9) ];
d2x = d2x(:, 10:end-9) + [zeros(length(X), N-9) d2x(:, 1:9)] + [d2x(:, end-8:end) zeros(length(X), N-9) ];
d3x = d3x(:, 10:end-9) + [zeros(length(X), N-9) d3x(:, 1:9)] + [d3x(:, end-8:end) zeros(length(X), N-9) ];

d0y = d0y(:, 10:end-9) + [zeros(length(X), N-9) d0y(:, 1:9)] + [d0y(:, end-8:end) zeros(length(X), N-9) ];
d1y = d1y(:, 10:end-9) + [zeros(length(X), N-9) d1y(:, 1:9)] + [d1y(:, end-8:end) zeros(length(X), N-9) ];
d2y = d2y(:, 10:end-9) + [zeros(length(X), N-9) d2y(:, 1:9)] + [d2y(:, end-8:end) zeros(length(X), N-9) ];
d3y = d3y(:, 10:end-9) + [zeros(length(X), N-9) d3y(:, 1:9)] + [d3y(:, end-8:end) zeros(length(X), N-9) ];

delta_functions = {d0x, d1x, d2x, d3x, d0y, d1y, d2y, d3y};
end
