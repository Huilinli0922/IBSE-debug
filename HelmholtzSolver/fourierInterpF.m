function y1 = fourierInterpF(pp,x1)
%	fourierInterpF  take spectral interpolation from pp = fft(u) to x1,y1
%
%       J. Mac Huang, 06/18/2016
%       final version 04/30/2023
%       AML DT, sculpt candy project
%
%   Inputs:  pp           - Fourier coefficients of interpolant
%            x1           - new grid
%
%   Outputs: y1           - interpolated value on x1
%

N = length(pp);

j1 = [0:N/2 -(N/2-1):-1]'*2*pi;

yt = pp;
y1 = 0*x1;
for k = 1:length(x1)
y1(k) = real(sum(yt.*exp(1i*j1*x1(k))))/N;
end

end