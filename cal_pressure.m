clc
clf
close all
clear all


server = 'office';    
workingFolder = '.';  
movieFolder = workingFolder;

filename = [workingFolder '/' server '_ice_' '03-Aug-2023' '.mat'];
% saveint = 1000;                         % how many steps until disk save
% saveskip = 100;                         % how many steps until memory save

path=['./',filename];
load(path);
pressure=[];
shear=[];
force=[];
position=[];
fx=[];
N=2;
surf_p = cell(1, N);  % Initialize a cell array
surf_tau = cell(1, N);
f_x = cell(1, N);
nx = cell(1, N);
ny = cell(1, N);
diag_stuff = cell(1, N);
x = cell(1, N);


for i = 1:N
    % data{i} = eval(['x', sprintf('%03d', i)]);
    surf_p{i} = eval(['surf_p', sprintf('%03d', i)]);
    surf_tau{i} = eval(['surf_tau', sprintf('%03d', i)]);
    f_x{i} = eval(['f_x', sprintf('%03d', i)]);
    nx{i} = eval(['nx', sprintf('%03d', i)]);
    ny{i} = eval(['ny', sprintf('%03d', i)]);
    diag_stuff{i} = eval(['diag_stuff', sprintf('%03d', i)]);
    x{i} = eval(['x', sprintf('%03d', i)]);

    
end


% surf_p={surf_p001,surf_p002};
% surf_tau={surf_tau001,surf_tau002};
% f_x={f_x001,f_x002};
% nx={nx001,nx002};
% ny={ny001,ny002};
% diag_stuff={diag_stuff001,diag_stuff002};
% x={x001,x002};

    % f_x = mean((-surf_p+diag_stuff).*nx+surf_tau.*ny);
    % f_y = mean((-surf_p-diag_stuff).*ny+surf_tau.*nx);

for i=1:N
for n=1:disk_save_int/save_skip


diag=diag_stuff{i}{n};
nxp=nx{i}{n};
tau=surf_tau{i}{n};
nyp=ny{i}{n};
p=surf_p{i}{n};
P=-mean(p.*nxp);
pressure=[pressure; P];
S=mean(diag.*nxp+tau.*nyp);
shear=[shear;S];
force=[force;P+S];

end
position=[position;x{i}(1:disk_save_int/save_skip)];
fx=[fx;f_x{i}(1:disk_save_int/save_skip)];
end

v_length=size(pressure,1);
time=(1:v_length)*dt*save_skip;

figure(1)
plot(time, force,'Color','r');
hold on
plot(time, fx(:,1),'Color','b');
hold on 
plot([0,400],[0,0],'color','k')

legend('force','fx')
title('total force')

figure(2)
plot(time, pressure,'Color','r')
hold on
plot(time, shear,'Color','b')
% ylim([-0.2,0.2])
legend('pressure','shear')
title('Comparison between pressure and shear')

figure(3)
plot(time, position)
title(' x position')

figure(4)
v=(position(2:end)-position(1:end-1))/(time(2)-time(1));
plot(time(1:end-1),v)
title ('velocity')
