function c = RK4(N, cold, u, v, dt)

C1 = - advective(N, u, v, cold); 
C2 = - advective(N, u, v, cold+1/3*dt*C1);
C3 = - advective(N, u, v, cold-1/3*dt*C1+dt*C2);
C4 = - advective(N, u, v, cold+dt*C1-dt*C2+dt*C3);

c = cold+dt/8*(C1+3*C2+3*C3+C4);


% C1 = - advective(N, u, v, cold); 
% C2 = - advective(N, u, v, cold+1/2*dt*C1);
% C3 = - advective(N, u, v, cold+1/2*dt*C2);
% C4 = - advective(N, u, v, cold+dt*C3);
% 
% c = cold+dt/6*(C1+2*C2+2*C3+C4);
end