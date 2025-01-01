function I = fourierInt_mf(u)
%	fourierInt_mf  take spectral integration of 1D periodic function that 
%                  mean free
%
%       J. Mac Huang, 06/18/2016
%       final version 07/03/2016
%       AML DT, sculpt candy project
%
%   Inputs:  u            - function
%
%   Outputs:
%            I            - int_0^s u ds
%

N = length(u);
if size(u,2)~=1
    u = u';
end
j1 = [0:N/2 -(N/2-1):-1]'*2*pi;
ubar = fft(u); 
ubar = ubar./(1j*j1); ubar(1) = 0;
I = real(ifft(ubar));
I = I - I(1);

end
