function v = Laplacian(u)
%	 Laplacian    take 2D laplacian of u
%
%       J. Mac Huang, 06/18/2016
%       final version 06/19/2016
%       AML DT, sculpt candy project
%
%   L equation: u_xx + u_yy - sigma*u = v
%

N = sqrt(length(u));
u = reshape(u,N,N);

j1 = repmat([0:N/2 -(N/2-1):-1]',1,N);
j2 = repmat([0:N/2 -(N/2-1):-1],N,1);
vbar = (- (j1.^2 + j2.^2)).*fft2(u);


v = real(ifft2(vbar));

v = v(:);
end