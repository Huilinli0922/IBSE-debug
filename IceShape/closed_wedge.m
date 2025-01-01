function [theta,X,Y,L]=closed_wedge(s,alpha,l0,smooth_parameter,height_ratio)
a=alpha;
b=pi/2-alpha;
L=(1+tan(b)+1/sin(a))*height_ratio;


% [X,Xn,theta] = computeBoundary(a,b,L,N); %compute boundary
N=length(s);
ds = s(2)-s(1);
h=smooth_parameter;
sigmoid = @(s0) 0.5*(1-tanh(-(s-s0)/h)); %smoothing function
% Solve for total arclength of each segment
options = optimset('TolFun',1e-15,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','none');
F_handle = @(val) F(val,a,b,sigmoid); %function to pass to fsolve
val = fsolve(F_handle,[1/4 1/4],options); %solve
[~,theta] = F(val,a,b,sigmoid); %compute tangent angle

% Fourier modes
if mod(N,2) == 0 %even
    iks_inv = 1./(2*pi*1i*[0:(N/2-1),0,(-N/2+1):-1]');
    iks_inv([1 N/2+1]) = 0;
elseif mod(N,2) == 1 %odd
    iks_inv = 1./(2*pi*1i*[0:(N-1)/2 -(N-1)/2:-1]');
    iks_inv(1) = 0;
end

% Solve for body
Xs = [cos(theta) sin(theta)]; %tangent vector
Xn = [Xs(:,2) -Xs(:,1)]; %normal vector
location = L*real(ifft(iks_inv.*fft(Xs))); %body
X=location(:,1);
Y=location(:,2);

% Check error
ex = sum(Xn(:,1))*ds;
ey = sum(Xn(:,2))*ds;
fprintf(" error int(nx) = %1.4e\n error int(ny) = %1.4e\n",ex,ey)

% if plot_initial==1
%     close all
% 
%     subplot(2,1,1)
%     plot(s,theta)
%     xlabel('$s$'),ylabel('$\theta$')
% 
%     subplot(2,1,2)
%     plot(X([1:end 1]),Y([1:end 1]),'o-')
%     hold on
%     quiver(X(:),Y(:),Xn(:,1),Xn(:,2))
%     axis equal
%     xlabel('$x$'),ylabel('$y$')
% 
%     title(['error \int nx=', num2str(ex), '; error \int nx=',num2str(ey) ])
% else
% end


end


function [err,theta] = F(val,a,b,sigmoid)

% Get arclength of each segment
L2 = val(1); %segment 2
L3 = val(2); %segment 3
L1 = 1/2-(L2+L3)/2; %segment 1 (same as segment 4)

% Cumulative arclength parameter
s1 = L1;
s2 = L1+L2;
s3 = L1+L2+L3;

% Evaluate tangent angle
theta = pi/2+(a-pi)*(sigmoid(s1)-sigmoid(s2)) + ...
    (-3*pi/2)*(sigmoid(s2)-sigmoid(s3)) + ...
    (-2*pi)*(sigmoid(s3));

% Compute error in endpoint
err = [mean(cos(theta));mean(sin(theta))];

end

