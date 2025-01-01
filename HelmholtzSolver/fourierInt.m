function [I,endVal] = fourierInt(u)
%	fourierInt  take spectral integration of 1D periodic function that is
%               not mean free
%
%       J. Mac Huang, 06/18/2016
%       final version 07/03/2016
%       AML DT, sculpt candy project
%
%   Inputs:  u            - function
%
%   Outputs:
%            I            - int_0^s u ds
%            endVal       - int_0^1 u ds
%
N = length(u);
if size(u,2)~=1
    u = u';
end
j1 = [0:N/2 -(N/2-1):-1]'*2*pi;
ubar = fft(u); 
meanu = mean(u);
ubar = ubar./(1j*j1); ubar(1) = 0;
I = real(ifft(ubar));
I = I - I(1) + meanu*(0:N-1)'/N;
endVal = meanu;

end
