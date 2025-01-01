% plot & print & save each step

% whether to print the solver messages
if ifprint
    if max(iter_c(2),iter_u(2))>iter_max
        renewflag2 = 1;
        disp( '  max iteration reached, SC will be renewed')
    end
    disp( ['    GMRES converged with ' num2str(max(iter_c(2),iter_u(2))) ...
        ' iterations'])
    disp( ['  timestep ' num2str(k) ', remaining length L/L_0= ' ...
        num2str(boundary_data{8}/L0*100, '%.2f') ' %'])
else
end

% plotting
if ifplot
    if mod(k, plotint) == 0
        % generate movie
        XOsharp = XO>0.5;

        figure(1)
        hh = pcolor(reshape(x,N,N),reshape(y,N,N),...
            reshape((c_wall-c)/(c_wall-c_ice).*(XOsharp)+XE,N,N));

        set(hh, 'EdgeColor','none')
        colormap(cmap)
        shading interp %smooth solution
        caxis([0,1]);
        colorbar
        hold on
        plot([X_np; X_np(1)],[Y; Y(1)], 'color', [0.3010 0.7450 0.9330], ...
            'lineWidth',2)
        plot([X_np; X_np(1)]-2*pi,[Y; Y(1)], 'color', [0.3010 0.7450 0.9330], ...
            'lineWidth',2)
        plot([X_np; X_np(1)]+2*pi,[Y; Y(1)], 'color', [0.3010 0.7450 0.9330], ...
            'lineWidth',2)
        hold on

        if plotpressure
            quiver(X_np,Y,pressure_field(:,1)/(2e6),pressure_field(:,2)/(2e6),...
                'AutoScale', 'off', 'color', 'r', 'lineWidth',1)
            hold on
        end


        if plotshear
            quiver(X_np(1:2:end),Y(1:2:end),shear_field(1:2:end,1),...
                shear_field(1:2:end,2),'AutoScale', 'off', 'color', 'g', 'lineWidth',1);
            hold on
        end

        %plot([0; 2*pi], 0*[0; 2*pi] + pi-pi/gamma, 'k', 'lineWidth',2 )
        %plot([0; 2*pi], 0*[0; 2*pi] + pi+pi/gamma, 'k', 'lineWidth',2 )
        axis equal
        %axis([0,2*pi, y_lower,y_upper])
        axis off
        title(['current time: ' num2str(tt*tau/60, '%.2f') ' min, f_x= '...
            num2str(f_x, '%.2e') ', L/L_0= ' ...
            num2str(boundary_data{8}/L0*100, '%.2f') ' %'], 'FontSize', font_size)
        drawnow
        writegif(gifname , pp_plot , 30, 1 )
        pp_plot = pp_plot+1;
        hold off
    end

    if mod(k,plotint*mkmovint) == 0 && ifmovie
        system(['LD_LIBRARY_PATH=""   && ffmpeg -y -i ' gifname ...
            ' -b:v 4M ' movieName]);
        system(['LD_LIBRARY_PATH=""   && chmod o+r ' movieName]);
    end


   
     % if mod(k,plotint*mkmovint) == 0 && ifmovie
     %        system(['ffmpeg -y -i ' gifname ...
     %            ' -crf 10 -b:v 1000K ' movieName]);
     %        system(['chmod o+r ' movieName]);
     % end

end

% save variables
if ifsave
    if mod(k, disk_save_int) ~= 0
        if mod(mod(k, disk_save_int), save_skip) == 0
            c_save{save_num} = c;
            XO_save{save_num} = XO;
            theta_save{save_num} = theta;
            L_save(save_num) = L;
            f_x_save(save_num)  = f_x;
            surf_p_save{save_num}  = surf_p;
            surf_tau_save{save_num}  = surf_tau;
            flow_save{save_num,1} = u;
            flow_save{save_num,2} = v;
            flow_save{save_num,3} = p;

            y0_save(save_num) = y0;
            tt_save(save_num) = tt;

            nx_save{save_num}  = nx;
            ny_save{save_num}  = ny;
            diag_stuff_save{save_num}  = diag_stuff;
            save_num = save_num+1;
        end
    else
        c_save{save_num} = c;
        XO_save{save_num} = XO;
        theta_save{save_num} = theta;
        L_save(save_num) = L;
        f_x_save(save_num)  = f_x;
        surf_p_save{save_num}  = surf_p;
        surf_tau_save{save_num}  = surf_tau;
        flow_save{save_num,1} = u;
        flow_save{save_num,2} = v;
        flow_save{save_num,3} = p;
        y0_save(save_num) = y0;
        tt_save(save_num) = tt;
        nx_save{save_num}  = nx;
        ny_save{save_num}  = ny;
        diag_stuff_save{save_num}  = diag_stuff;
        %             bd_velocity_save{save_num}=bd_velocity; %[u_bd,v_bd,Vx,Vy]
        disp('    saving...')

        %create variable names
        cvar = ['c',sprintf('%03d',pp_save)];
        XOvar = ['XO',sprintf('%03d',pp_save)];
        Lvar = ['L',sprintf('%03d',pp_save)];
        fx_var = ['f_x',sprintf('%03d',pp_save)];
        surf_p_var = ['surf_p',sprintf('%03d',pp_save)];
        surf_tau_var = ['surf_tau',sprintf('%03d',pp_save)];
        thetavar = ['theta',sprintf('%03d',pp_save)];
        flowvar = ['flow',sprintf('%03d',pp_save)];
        y0var = ['y',sprintf('%03d',pp_save)];
        ttvar = ['tt',sprintf('%03d',pp_save)];
        nxvar = ['nx',sprintf('%03d',pp_save)];
        nyvar = ['ny',sprintf('%03d',pp_save)];
        diag_stuffvar=['diag_stuff',sprintf('%03d',pp_save)];
        % bd_velocityvar = ['bd_velocity',sprintf('%03d',pp_save)];

        eval([cvar, '= c_save;']);
        eval([XOvar, '= XO_save;']);
        eval([thetavar, '= theta_save;']);
        eval([flowvar, '= flow_save;']);
        eval([Lvar, '= L_save;']);
        eval([fx_var, '= f_x_save;']);
        eval([surf_p_var, '= surf_p_save;']);
        eval([surf_tau_var, '= surf_tau_save;']);
        eval([y0var, '= y0_save;']);
        eval([ttvar, '= tt_save;']);
        eval([nxvar, '= nx_save;']);
        eval([nyvar, '= ny_save;']);
        eval([diag_stuffvar, '= diag_stuff_save;']);
        save(filename, cvar, thetavar, Lvar, fx_var, surf_p_var, ...
            surf_tau_var, flowvar , XOvar, y0var, ttvar,nxvar,nyvar,...
            diag_stuffvar, '-append');
        pp_save = pp_save+1;
        save_num = 1;
    end
end

if ifsave_dt
    if mod(k, disk_save_int2) ~= 0
        total_p_save(save_num2,:) = pressure;
        total_shear_save(save_num2,:) = shear;
        x0_save(save_num2) = x0;
        u_ice_save(save_num2) = u_ice;
        s_save(save_num2) = S0;
        save_num2 = save_num2+1;
    else
        total_p_save(save_num2,:) = pressure;
        total_shear_save(save_num2,:) = shear;
        x0_save(save_num2) = x0;
        u_ice_save(save_num2) = u_ice;
        s_save(save_num2) = S0;

        %create variable names
        total_pvar = ['total_p',sprintf('%03d',pp_save2)];
        total_shearvar = ['total_shear',sprintf('%03d',pp_save2)];
        x0var = ['x',sprintf('%03d',pp_save2)];
        u_icevar = ['u_ice',sprintf('%03d',pp_save2)];
        svar = ['s',sprintf('%03d',pp_save2)];

        eval([total_pvar, '= total_p_save;']);
        eval([total_shearvar, '= total_shear_save;']);
        eval([x0var, '= x0_save;']);
        eval([u_icevar, '= u_ice_save;']);
        eval([svar, '= s_save;']);

        save(filename, total_pvar, total_shearvar,x0var,u_icevar,svar, '-append');
        pp_save2 = pp_save2+1;
        save_num2 = 1;
    end
end
