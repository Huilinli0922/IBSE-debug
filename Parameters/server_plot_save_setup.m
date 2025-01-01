% server&plot set up

%% server location
if location == 0
    server = 'cims';
    cores = 30;
    workingFolder = 'data/';
    ifvisible = 0;
elseif location == 1
    server = 'office';
    cores = 48;
    workingFolder = '~/Scratch/iceboat';
    ifvisible = 1;
elseif location == 2
    server = 'tiffany';
    cores = 3;
    workingFolder = 'data/';
    ifvisible = 1;
elseif location == 3
    server = 'mac';
    cores = 4;
    workingFolder = 'D:\machuang\Work\iceboat';
    ifvisible = 1;
elseif location == 4
    server = 'nyush_hpc';
    cores = 32;
    workingFolder = '~/Scratch/iceboat';
    ifvisible = 0;
elseif location == 5
    server = ['nyush_deg' num2str(round(alpha*180/pi))];
    cores = 32;
    workingFolder = 'data/';
    ifvisible = 0;
end

movieFolder = workingFolder;            % where to save movie
warning('off','all')                    % ignore GMRES warnings

%% Solver work path, plotting and saving behaviors
ifsave = 1;                             % whether to save
ifsave_dt = 1;                          % whether to save every other step
filename = [workingFolder '/' date '_' server  '.mat'];


ifprint =1;                             % whether to print
ifplot = 1;                             % whether to plot
ifmovie = 0;                            % whether to generate mp4 movie
%     Note: requires FFMPEG
plotint = 10;                           % skip of plot
mkmovint =30;                         % skip of making movie
hidpi = 0;                              % whether the display is hidpi
%     Note: set to 0 on server

plotpressure=1;
plotshear=0;

% path for gif and mp4 movies
gifname = [workingFolder '/' date '_' server '.gif'];
movieName = [movieFolder '/' date '_' server '.mp4'];

% set colormap
cmap = brewermap(256, 'blues');

% prepare for saving and plotting
if ifplot
    Fig = figure(1);
    if ~ifvisible
        set(Fig,'visible','off');
    end
    set(Fig, 'Position', [1 1 800 800])
    if hidpi
        font_size = 16*settings().matlab.desktop.DisplayScaleFactor.PersonalValue;
    else
        font_size = 16;
    end
end

if ifsave
    save(filename);

    theta_save = cell( floor(disk_save_int/save_skip+1),1);
    L_save = zeros( floor(disk_save_int/save_skip+1),1);
    f_x_save = L_save; f_tau_save = L_save;
    c_save = cell( floor(disk_save_int/save_skip+1),1);
    XO_save = cell( floor(disk_save_int/save_skip+1),1);
    flow_save = cell( floor(disk_save_int/save_skip+1),3);
    y0_save = zeros( floor(disk_save_int/save_skip+1),1);
    tt_save = zeros( floor(disk_save_int/save_skip+1),1);
    surf_p_save = cell( floor(disk_save_int/save_skip+1),1);
    surf_tau_save = cell( floor(disk_save_int/save_skip+1),1);
    nx_save = cell( floor(disk_save_int/save_skip+1),1);
    ny_save = cell( floor(disk_save_int/save_skip+1),1);
    diag_stuff_save = cell( floor(disk_save_int/save_skip+1),1);
    % variables saved every other step
    total_p_save =zeros(disk_save_int2,2);
    total_shear_save =zeros(disk_save_int2,2);
    x0_save = zeros( disk_save_int2,1);
    u_ice_save= zeros( disk_save_int2,1);
    s_save=zeros( disk_save_int2,1);

else
end


%% Preparation
disp('    starting...')

% start parallel pool
pll = gcp('nocreate'); % If no pool, do not create new one.
if isempty(pll)
    parobj = parpool('local',cores);
end
%pctRunOnAll warning off


