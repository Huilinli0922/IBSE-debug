function imtheta = implicitTheta(theta, epsilon, L, dt,sigma)
nbdy0 = length(theta);
j1 = [0:nbdy0/2 -(nbdy0/2-1):-1]';

imtheta = real(ifft(fft(theta)./(1 + epsilon*dt*sigma/L^2*j1.^2)));
end
