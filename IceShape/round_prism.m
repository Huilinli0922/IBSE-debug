
function [theta,X,Y,L] = round_prism(s, alpha, h, r );
% %%smooth boundary
% clear all
% close all
% format long


% L0=1/sin(alpha)+1+1/tan(alpha);

%   L           -total arc lenghth
%   X,Y         -position of boundary 
%   theta       -tangent angle with y^ direction

nbdy=length(s);

l1=(h-2*r)*cot(alpha);
l2=(pi-alpha)*r;
l3=(h-2*r)/sin(alpha);
l4=(pi/2+alpha)*r;
l5=h-2*r;
l6=pi/2*r;
l=[l1,l2,l3,l4,l5,l6];
L=sum(l);

ds=sum(l)/nbdy;
dtheta=ds/r;
n=zeros(nbdy,2); X=zeros(nbdy,1); Y=zeros(nbdy,1); theta=zeros(nbdy,1);

N1=floor(l(1)/ds);
N2=floor(sum(l(1:2))/ds);
N3=floor(sum(l(1:3))/ds);
N4=floor(sum(l(1:4))/ds);
N5=floor(sum(l(1:5))/ds);
N6=floor(sum(l(1:6))/ds);

for i=1:N1
    theta(i)=0;
    n(i,:)=[0,1];
    X(i)=ds*i;
    Y(i)=r;
end

for i=N1+1
    b1=l1-N1*ds;
    b2=ds-b1;
    t1=b2/r;
    theta(i)=t1;
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    X(i)=l1+sin(theta(i))*r;
    Y(i)=cos(theta(i))*r;
end

for i=(N1+2):N2

    theta(i)=t1+dtheta*(i-N1-1);
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    X(i)=l1+sin(theta(i))*r;
    Y(i)=cos(theta(i))*r;
end

for i=N2+1
    theta(i)=pi-alpha;
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    t3=pi-alpha-theta(N2);
    b3=t3*r;
    b4=ds-b3;
    Bp=[l1+sin(theta(i))*r,cos(theta(i))*r];
    X(i)=Bp(1)+b4*sin(theta(i)+pi/2);
    Y(i)=Bp(2)+b4*cos(theta(i)+pi/2);
end

for i=N2+2:N3
    theta(i)=pi-alpha;
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    X(i)=Bp(1)+(b4+(i-N2-1)*ds)*sin(theta(i)+pi/2);
    Y(i)=Bp(2)+(b4+(i-N2-1)*ds)*cos(theta(i)+pi/2);
end

for i=N3+1
    b5=l1+l2+l3-ds*N3;
    b6=ds-b5;
    t6=b6/r;
    theta(i)=t6+pi-alpha;
    C=[0,-l5];
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    X(i)=C(1)+r*sin(theta(i));
    Y(i)=C(2)+r*cos(theta(i));
end

for i=(N3+2):N4
    theta(i)=t6+pi-alpha+dtheta*(i-N3-1);
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    X(i)=C(1)+r*sin(theta(i));
    Y(i)=C(2)+r*cos(theta(i));
end

for i=N4+1
    t7=3*pi/2-theta(N4);
    b7=t7*r;
    b8=ds-b7;
    theta(i)=3/2*pi;
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    Cp=C+r*[sin(3*pi/2),cos(3*pi/2)];
    X(i)=Cp(1)+b8*sin(2*pi);
    Y(i)=Cp(2)+b8*cos(2*pi);
end

for i=N4+2:N5
    theta(i)=3*pi/2;
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    % n(i,:)=[-1,0];
    X(i)=Cp(1)+(b8+ds*(i-N4-1))*sin(2*pi);
    Y(i)=Cp(2)+(b8+ds*(i-N4-1))*cos(2*pi);
end

for i=N5+1
    b9=sum(l(1:5))-N5*ds;
    b10=ds-b9;
    t10=b10/r;
    theta(i)=t10+3*pi/2;
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    A=[0,0];
    X(i)=A(1)+r*sin(theta(i));
    Y(i)=A(2)+r*cos(theta(i));
end

for i=N5+2:N6
    theta(i)=t10+3*pi/2+dtheta*(i-N5-1);
    n(i,:)=[sin(theta(i)),cos(theta(i))];
    X(i)=A(1)+r*sin(theta(i));
    Y(i)=A(2)+r*cos(theta(i));
end
theta=-theta+pi/2;

X=X-X(1);
Y=Y-Y(1);
% plot(X,Y,'linewidth',2)
% hold on 

% 
% quiver(X(:),Y(:),n(:,1),n(:,2));
% axis equal
% legend('original','forier integral','normal vector')
% xlim([-1,2])
% ylim([-1.5,0.5])
% 
% disp(['int nx=',num2str(mean(n(:,1)))])
% disp(['int ny=',num2str(mean(n(:,2)))])
% disp(['int nx*Y=',num2str(mean(n(:,1).*Y))])
end


% arc=sqrt((Xs(2:end)-Xs(1:(end-1))).^2+(Ys(2:end)-Ys(1:end-1)).^2);
% arc_sum=sum(arc);
% 
% arc2=sqrt((X(2:end)-X(1:(end-1))).^2+(Y(2:end)-Y(1:end-1)).^2);
% arc2_sum=sum(arc2);

% function I = fourierInt_mf2(u)
% N = length(u);
% 
% if size(u,2)~=1
%     u = u';
% end
% 
% j1 = [0:N/2 -(N/2-1):-1]'*2*pi;
% ubar = fft(u); 
% init=ubar(1);
% ubar = ubar./(1j*j1); ubar(1) = 0;
% % ubar(1)=init;
% I = real(ifft(ubar));
% I = I - I(1);
% 
% 
% 
% end


% function I = fourierInt2(u)
% N = length(u);
% if size(u,2)~=1
%     u = u';
% end
% j1 = [0:N/2 -(N/2-1):-1]'*2*pi;
% ubar = fft(u); 
% meanu = mean(u);
% % meanu = ubar(1);
% ubar = ubar./(1j*j1); ubar(1) = 0;
% I = real(ifft(ubar));
% I = I - I(1) + meanu*(0:N-1)'/N;
% endVal = meanu;
% end
