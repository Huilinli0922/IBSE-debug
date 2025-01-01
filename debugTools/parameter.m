%calculating parameters used for ice melting problem 
% Tiffany Li at 06/25/2023
% all the units are interational standard units

T=20; %temperature-K
nu=1e-6; %kinematic viscosity of water- m^2/s
mu=1e-3; %dynamic viscosity Ns/m^2
a_v=207e-6;%thermal expansion coeff
g=9.8;
L=0.07;
T_m=0; %melting temperature in celcius degree
T_infty=20; %ambient temperature in celcius degree
dT=20; %temperature difference between ice and water
K=0.143e-6; %thermal diffusivity
k=0.598; %thermal conductivity 
c_p=4184; %specific heat
l=334e3; % latent heat 

Gr=L^3*a_v*dT*g/(nu^2); %Grashof number
Re=sqrt(Gr); %reynolds number 
Pr=c_p*mu/k; %prandtl number
St=c_p*dT/l; %Stephan number 
Pe=Re*Pr; %Peclet number 
Ra=Re*Pe; %Rayleigh number
