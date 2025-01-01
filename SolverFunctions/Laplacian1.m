function v = Laplacian1(u)
%	 Laplacian1    take 1D laplacian of u
%
%       J. Mac Huang, 06/18/2016
%       final version 06/19/2016
%       AML DT, sculpt candy project
%
%

N = (length(u));

j1 =[0:N/2 -(N/2-1):-1]';
vbar = (- (j1.^2)).*fft(u - mean(u));

v = real(ifft(vbar));

v = v(:);
end