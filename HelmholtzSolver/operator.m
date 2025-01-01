function v = operator(u, OP)
%	operator compute v = OP*u with Fourier spectral method, here OP
%           is an elliptic operator.
%
%       J. Mac Huang, 06/18/2016
%       final version 06/19/2016
%       AML DT, sculpt candy project
%
N = sqrt(length(u));
v = real(ifft2(fft2(reshape(u,N,N)).*OP));
v = v(:);
end
