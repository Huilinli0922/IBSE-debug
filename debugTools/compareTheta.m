function err = compareTheta(theta1,theta2)
N1 = length(theta1);
N2 = length(theta2);
s1 = (0:N1-1)'/N1;
s2 = (0:N2-1)'/N2;

if N1 < N2
    diff = theta1-interp1(s2,theta2,s1,'pchip');
    err = norm(diff,inf)/norm(theta2,inf);
else
    diff = theta2-interp1(s1,theta1,s2,'pchip');
    err = norm(diff,inf)/norm(theta1,inf);
end
end
