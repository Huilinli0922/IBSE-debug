function R0 = getR0(lambda)

y = @(x) (x.^2/4.*expint(x.^2/4).*exp(x.^2/4)-lambda);


R0 = fsolve(y,lambda);
end