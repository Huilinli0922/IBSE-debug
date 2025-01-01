function saveStefan(k, c, theta, L, dt0, N)
eval(['c' num2str(k) '=c;']);
eval(['theta' num2str(k) '=theta;']);
eval(['L' num2str(k) '=L;']);
eval(['dt' num2str(k) '=dt0;']);
eval(['N' num2str(k) '=N;']);
save(['C' num2str(k)],['c' num2str(k)], ['theta' num2str(k)],['L' num2str(k)], ['dt' num2str(k)] )
end