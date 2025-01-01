function [f, y_lower, y_upper] = damping_mask(y, y_lower, y_upper, width)
%   damping_mask   define the damping area for the absorbing boundary cond. 
%
%       J. Mac Huang, 07/17/2023
%       AML DT, iceboat project
%
%   Inputs: boundary_data      - boundary data generated by boundary_info
%           damping_strength   - strength of the damping function
%           sigma              - width of the damping region
%           filter_loc         - location of filter on lower bound (< pi)
%
%   Output: f                  - damping function
%           Yw1_new            - location of the fictitious wall at buttom
%           Yw2_new            - location of the fictitious wall at top
%

if y_upper + width/2 > 2*pi
    y_upper = 2*pi - width/2;
end

if y_lower - width/2 < 0
    y_lower = width/2;
end

h_upper = smooth_step_func(y - y_upper, width);
h_lower = smooth_step_func(y_lower - y, width);

f = h_lower+h_upper;

end