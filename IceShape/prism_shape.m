function [theta, L, candy_side] = prism_shape(s, alpha, l0, smooth_parameter)
%       prism_shape   define the geometry of a ice prism with right angle
%       and fixed cathetus
%  
%
%       J. Mac Huang, 09/27/2022
%       AML DT, candy boat project
%       Tiffany Li, modified at 06/22/2023
%
%   Inputs:  s                  - normalized arclength in [0,1)
%            alpha              - candy boat angle
%            l0                 - length of the candy
%            smooth_parameter   - parameter to round the corner, smaller
%                              value means sharper corners
%
%   Outputs: theta              - tangent angle of the boat surface
%            L                  - total arclength of the boat
%            candy_side         - a smooth characteristic function
%                            indicating the area of boat covered by candy
%

% arclength for 2 corners of the candy boat
s1 = 1/(cos(alpha)+sin(alpha)+1)*cos(alpha);
s2 = 1/(cos(alpha)+sin(alpha)+1)*(cos(alpha)+1);


% define the triangle of the boat shape
theta = pi/2*(s<s1) + (alpha-pi/2)*(s<s2).*(s>=s1) ...
                    - pi.*(s>=s2);
% theta(end/2+1) = (theta(end/2)+theta(end/2+2))/2;

% compute total arclength
L = (1+cos(alpha)+sin(alpha))*l0;

% smooth the profile of theta through Gaussian filter
theta = gaussianF( theta+s*2*pi, smooth_parameter ) - s*2*pi;

% define the smooth characteristic function for the candy surface, Wendland
% fundtion of C2 is used here.
candy_side = smooth_step_func(s-s1,smooth_parameter/2).*...
             smooth_step_func(s2-s,smooth_parameter/2);
end
        
  
  