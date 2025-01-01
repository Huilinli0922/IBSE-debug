% create all the save variables

theta_save = cell( floor(disk_save_int/save_skip+1),1);
L_save = zeros( floor(disk_save_int/save_skip+1),1);
f_x_save = L_save; f_tau_save = L_save;
c_save = cell( floor(disk_save_int/save_skip+1),1);
XO_save = cell( floor(disk_save_int/save_skip+1),1);
flow_save = cell( floor(disk_save_int/save_skip+1),3);
x0_save = zeros( floor(disk_save_int/save_skip+1),1);
y0_save = zeros( floor(disk_save_int/save_skip+1),1);
tt_save = zeros( floor(disk_save_int/save_skip+1),1);
surf_p_save = cell( floor(disk_save_int/save_skip+1),1);
surf_tau_save = cell( floor(disk_save_int/save_skip+1),1);
u_ice_save= zeros( floor(disk_save_int/save_skip+1),1);
nx_save = cell( floor(disk_save_int/save_skip+1),1);
ny_save = cell( floor(disk_save_int/save_skip+1),1);
diag_stuff_save = cell( floor(disk_save_int/save_skip+1),1);
% bd_velocity_save=cell( floor(disk_save_int/save_skip+1),1); %bd_velocity={u_bd,v_bd,Vx,Vy};
% variables saved every other step
total_p_save =zeros(disk_save_int2,2);
total_shear_save =zeros(disk_save_int2,2);