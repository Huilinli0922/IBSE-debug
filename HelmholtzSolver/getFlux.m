function flux = getFlux(normD,boundary_size, c)
%	getFlux  get flux from concentration field
%
%       J. Mac Huang, 06/18/2016
%       final version 07/03/2016
%       AML DT, sculpt candy project
%
%   Inputs:  normD          - boundary gradient operator
%            boundary_size  - size of the boundary
%            c              - concentration field
%
%   Outputs: flux           - flux at boundary
%

flux = normD*c; flux = flux(1:boundary_size);


end
