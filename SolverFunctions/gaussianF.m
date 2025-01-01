function dv = gaussianF( v, sigma )
% Gaussian filter
% INPUTS: vector of data,  Gaussian strength
%       Nick Moore
%
threshould = 1e-10;

N = length(v);
if( mod(N,2)==0 )
    k = ((-N/2):(N/2-1))';
else
    k = (-(N-1)/2:(N-1)/2)'; 
end

gv = exp(-0.5*abs(k).^2*sigma^2);
fv = fftshift(fft(v-mean(v))); fv(abs(fv)<threshould)=[0.];
fv = fv.*gv;
dv = real(ifft(ifftshift(fv)))+mean(v);
end

