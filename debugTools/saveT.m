function saveT(k,c, XO,omega,theta,L, dt0)
eval(['c' num2str(k) '=c;']);
eval(['omega' num2str(k) '=omega;']);
eval(['theta' num2str(k) '=theta;']);
eval(['L' num2str(k) '=L;']);
eval(['XO' num2str(k) '=XO;']);
eval(['dt' num2str(k) '=dt0;']);
save(['T' num2str(k)], ['XO' num2str(k)],['c' num2str(k)], ['omega' num2str(k)], ['theta' num2str(k)],['L' num2str(k)], ['dt' num2str(k)] )
end