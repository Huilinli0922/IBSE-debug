function boundary_data = boundaryinfo(theta, L, x0, y0, x, y, gamma)
%	boundaryinfo     get boundary information
%
%       J. Mac Huang, 06/18/2016
%       final version 05/20/2016
%       AML DT, sculpt candy project
%
%
%   Inputs:  theta        - tangent angle on the boundary Gamma
%            L            - total arclength of Gamma
%            x0,y0        - s=0 coordinate of the boundary Gamma
%            x, y         - meshgrid of the whole computational domain
%            gamma        - the aspect ratio of the solving regian
%                           (width/height)
%            nwall        - number of nodes on the wall, shold be
%                           approximately N/2
%
%   Outputs: X, Y         - boundary nodes on gamma
%            nx, ny       - outward normal direction of the boundary
%            dS           - metric on the boundary, dS = || dX/ds||
%            XE           - exterior
%            XEw          - exterior of the wall (so interior XO = ~XE.*(~XEw))
%            Xw1, Yw1     - nodes on the lower wall
%            Xw2, Yw2     - nodes on the upper wall
%

Nx = sqrt(length(x));    % total number of grid pts along x
n = 2;                   % half number of grids in the smooth transition

% get the rescaled arclength s and metric dS.
N = length(theta); dx = 2*pi/(sqrt(length(x)));
s = (0:N-1)'/N; dS = L/N+0*s;

% extend domain to accomodate for periodic boundary conditions
[x_e,y_e] = meshgrid(((-Nx+1):(2*Nx))/Nx*2*pi, (1:Nx)/Nx*2*pi);
x_e = x_e(:); y_e = y_e(:);

x0_np = x0;         % keep a record of nonperiodic value of x0
x0 = mod(x0, 2*pi);
X = x0 + L*fourierInt_mf(sin(theta)); X = X - (X(1)-x0);
Y = y0 - L*fourierInt_mf(cos(theta)); Y = Y - (Y(1)-y0);

% define the outward normal direction of boundary Gamma
nx = cos(theta);
ny = sin(theta);

% wave numbers on Gamma
j1 = [0:N/2-1 0 -(N/2-1):-1]'*2*pi;


% Fourier transformation
theta_hat = fft(theta + 2*pi*s);
theta_s_hat = 1i*j1.*fft(theta + 2*pi*s);
X_hat = fft(X);
Y_hat = fft(Y);
nx_hat = fft(nx);
ny_hat = fft(ny);

% preallocate XE
XE = x_e*0;

parfor j = 1:length(x_e)
    [dist, ind] = min((x_e(j) - X).^2 + (y_e(j) - Y).^2);
    dist = sqrt(dist);    % the distance between current grid to Gamma
    
    if dist/dx>n+2
        % if distance is far, no need for smooth
        % checking whether the grid is in or out of XE
        if nx(ind)*(x_e(j) - X(ind))+ny(ind)*(y_e(j) - Y(ind))<0
            XE(j) = 1;    % inside
        else
            XE(j) = 0;    % outside
        end
    else
        % if distance is near, it has to be computed through Newton method
        % s_min is the s location that minimize dist
        s_min = s(ind);         % initial guess
        res = 1; count = 1;     % initiate iteration
        X_min = 0; Y_min = 0;
        while abs(res) > 1e-8
            % evaluate theta, theta_s at s_min
            theta_min = fourierInterpF(theta_hat,s_min) - 2*pi*s_min;
            theta_s_min = fourierInterpF(theta_s_hat,s_min) - 2*pi;
            
            % evaluate X and Y at s_min
            X_min = fourierInterpF(X_hat,s_min);
            Y_min = fourierInterpF(Y_hat,s_min);
            
            % we want dist_s = 0 to find the minimal dist
            % first s derivative of dist
            dist_s = -L*cos(theta_min)*(Y_min-y_e(j)) + L*sin(theta_min)*(X_min-x_e(j));
            % second s derivative of dist
            dist_ss = L^2 + L*cos(theta_min)*(X_min-x_e(j))*theta_s_min + ...
                L*sin(theta_min)*(Y_min-y_e(j))*theta_s_min;
            
            % do Newton iteraction
            res = dist_s/dist_ss;
            s_min = s_min - res;
            
            count = count+1;
            if count>200
                warning('Newton method in boundary info did not converge!');
                break;
            end
        end
        
        % updated distance between (x(j), y(j)) and Gamma
        dist = sqrt((x_e(j)-X_min)^2 + (y_e(j)-Y_min)^2);
        
        nx_min = fourierInterpF(nx_hat,s_min);
        ny_min = fourierInterpF(ny_hat,s_min);
        side = - sign(nx_min*(x_e(j) - X_min)+ny_min*(y_e(j) - Y_min));
        % side is 1 in XE, 0 otherwise
        
        % uncomment here if you want to smooth the XE
        % XE(j) = wendland_2(n, dist/dx*side);
        
        
         % uncomment here if you don't want to smooth the XE
        if side >0
           XE(j) = 1;
        else
          XE(j) = 0;
        end
    end
end

% bring XE back from the extended periodic domain
XE = reshape(XE, Nx, 3*Nx);
XE = XE(:,1:Nx) + XE(:,Nx+1:2*Nx) + XE(:,2*Nx+1:3*Nx); XE = XE(:);
%XE = XE > 0.5;

% computhe XE for the wall, this is way easier as the wall shape is
% simple...

% comment here if you want to smooth the XEw
% XEw = wendland_2(n, (y-(pi+pi/gamma))/dx) + wendland_2(n, (-y+(pi-pi/gamma))/dx) ;

% uncomment here if you don't want to smooth the XEw

% indicator for physical domain
XO = (1-XE);

X_np = X;           % keep a nonperiodic X for plotting

% area of the ice
S0 = polyarea(interpft(X,10000),interpft(Y,10000));

X = mod(X, 2*pi);



boundary_data = {X, Y, nx, ny, XE, XO, ...
    theta, L, L, x0_np, y0, x, y, gamma, X_np, S0};
end
