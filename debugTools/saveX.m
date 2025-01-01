function saveX(k, c, XO, omega, psi, theta, L, dt0)
eval(['c' num2str(k) '=c;']);
eval(['omega' num2str(k) '=omega;']);
eval(['theta' num2str(k) '=theta;']);
eval(['L' num2str(k) '=L;']);
eval(['psi' num2str(k) '=psi;']);

eval(['XO' num2str(k) '=XO;']);
eval(['dt' num2str(k) '=dt0;']);
save(['X' num2str(k)], ['XO' num2str(k)],['c' num2str(k)],['psi' num2str(k)], ['omega' num2str(k)], ['theta' num2str(k)],['L' num2str(k)], ['dt' num2str(k)] )
end