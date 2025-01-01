function h = wendland_2(N,d)
%	wendland_2       generate Wenland function - regularized Heaviside
%                    function
%
%       J. Mac Huang, 06/18/2016
%       final version 06/19/2016
%       AML DT, sculpt candy project
%
%
%   Inputs:  N            - half number of grid points in the smooth region
%            d            - number of grid away from the center

%
%   Outputs: h            - regularized Heaviside function
%
h = 0*d;
for i = 1:length(d)
    x = abs(d(i))/N;
    if d(i)>0
        if d(i)>N
            h(i) = 1;
        else
            h(i) = 0.5+(x-28/9*x.^3+14*x.^5-224/9*x.^6+20*x.^7-8*x.^8+35/27*x.^9)*27/16;
        end
    else
        if d(i)<-N
            h(i) = 0;
        else
            h(i) = 0.5-(x-28/9*x.^3+14*x.^5-224/9*x.^6+20*x.^7-8*x.^8+35/27*x.^9)*27/16;
        end
    end
end
end