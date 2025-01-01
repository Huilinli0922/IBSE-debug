function err = compareSol(c1,XO1,c2,XO2)
C1 = c1.*(XO1==1); C2 = c2.*(XO2==1);
N1 = sqrt(length(C1));
N2 = sqrt(length(C2));

[X1, Y1] = meshgrid((1:N1)/N1*2*pi,(1:N1)/N1*2*pi);
[X2, Y2] = meshgrid((1:N2)/N2*2*pi,(1:N2)/N2*2*pi);


C1 = reshape(C1,sqrt(length(C1)),sqrt(length(C1)));
C2 = reshape(C2,sqrt(length(C2)),sqrt(length(C2)));
XO1 = reshape(XO1,N1,N1);
XO2 = reshape(XO2,N2,N2);

if N1 < N2
    diff = (C1-interp2(X2,Y2,C2,X1,Y1,'cubic')).*(XO1==1).*(interp2(X2,Y2,XO2,X1,Y1,'cubic')==1);
    diff = diff(:);
    err = norm(diff,inf)/norm(C2.*XO2, inf);
else
    diff = (C2-interp2(X1,Y1,C1,X2,Y2,'cubic')).*(XO2==1).*(interp2(X1,Y1,XO1,X2,Y2,'cubic')==1);
    diff = diff(:);
    err = norm(diff,inf)/norm(C1.*XO1, inf);
end


end
