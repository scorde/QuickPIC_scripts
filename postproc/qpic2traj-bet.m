%
% Matlab script to get electron trajectories 
% from QuickPIC output data and to compute betatron radiation
% 
% S. Corde, Jan 26, 2013
%
% Last update: S. Corde, Feb 16, 2013
%
% Only for QuickPIC binary qpic.e.0125
%


%% Defining input parameters

% update working dir to your QUICKPICSIM folder
working_dir = '~/Dropbox/SeB/Codes/sources/QuickPIC_scripts/';

% update data dir to your quickpic output folder (folder which contains the rp input file)
sim_number = 102;
if sim_number == 0
    datadir = '~/Dropbox/SeB/PostDoc/Projects/Electron-symmetrisation/QuickPIC_sim/qpic.e.0125_success/';
else
    datadir = ['~/Dropbox/SeB/PostDoc/Projects/Electron-symmetrisation/QuickPIC_sim/Sim_' num2str(sim_number) '/'];
end
% datadir = '~/Dropbox/SeB/PostDoc/Projects/Electron-symmetrisation/QuickPIC_sim/f/';

% Define a memory size to use for storing trajectories, in GB
RAM = 2;

% Indicate if tags are available, and if phi/psi potentials are available
tags_on = 0;
potential_on = 0;

% Parameters needed for betatron computation
home = '/Users/scorde/';
bet_executable = 'bet';
n_process = 8;

% start plot 3D timestep, -1: start from start
n_3D_timestep_start = -1;

% -1: go all the way to end
% n_3D_timestep_end = -1;
n_3D_timestep_end = 4240;



addpath(working_dir);

my_SI_params;
custom_cmap;

myfile_rpinput = [datadir 'rpinput'];
n0 = my_get_quickpic_param(myfile_rpinput, 'Plasma_Density');
TEND = my_get_quickpic_param(myfile_rpinput, 'TEND');
DT = my_get_quickpic_param(myfile_rpinput, 'DT');
DFPHA_BEAM = my_get_quickpic_param(myfile_rpinput, 'DFPHA_BEAM');

omega_p = 5.64e4*sqrt(n0);  % Plasma frequency in s-1

if(n_3D_timestep_start == -1)
  n_3D_timestep_start = 0;
end
if(n_3D_timestep_end == -1)
  n_3D_timestep_end = floor(TEND/DT);
end


npt = floor((n_3D_timestep_end-n_3D_timestep_start)/DFPHA_BEAM) + 1;
T = DFPHA_BEAM*DT*(npt-1)/omega_p;
npart = size(my_get_quickpic_phasespace(datadir, '01', sprintf('%.4d', n_3D_timestep_start+(npt-1)*DFPHA_BEAM)), 1);


disp(['Time step range is from ' num2str(n_3D_timestep_start, '%.4d') ' to ' num2str(n_3D_timestep_start+(npt-1)*DFPHA_BEAM, '%.4d')]);



%% Reading of the distribution files and writing of the trajectory file qp_traj

if npart*npt*6> RAM*1e9/8    
    npart_per_step = floor( RAM*1e9/8/(6*npt) ); % Define a number of particle to save no more than RAM GB of data at a time
    n_saving_step = floor(npart/npart_per_step);
    reste = npart - n_saving_step * npart_per_step;
    npart_save = 0;
    fid = fopen([datadir 'qp_traj'], 'wb');
    for i=1:n_saving_step
        istart = 1 + (i-1)*npart_per_step;
        istop = i*npart_per_step;
        npart_save = npart_save + save_qp_traj(istart, istop, n_3D_timestep_start, npt, datadir, fid);
    end
    npart_save = npart_save + save_qp_traj(npart-reste+1, npart, n_3D_timestep_start, npt, datadir, fid);
    fclose(fid);
else
    fid = fopen([datadir 'qp_traj'], 'wb');
    npart_save = save_qp_traj(1, npart, n_3D_timestep_start, npt, datadir, fid);
    fclose(fid);
end

Q = 3e3*npart_save/npart;

disp('Total number of macro-particles saved:');
disp(npart_save);
disp('Number of time steps saved:');
disp(npt);
disp('Total time T to input in the betatron code');
disp(T);
disp('Beam charge Q to input in the betatron code');
disp(Q);



%% Write the input file qp_bet_param.txt for the betatron code

distance_plasma_lanex = 23.27; % in m

param_header = sprintf('input_type = 1\nsync_limit = 1\nnpt_rec = 1\nnps = 0\nnps_spacing = 0\ntheta_X_max = 2.5\ntheta_Y_max = 2.5\nangular_resolution = 0.025\nsave_traj = 0\nsave_txt_d2W = 0\n');
param_file = [param_header, sprintf('traj_file = %s%sqp_traj\n', home, datadir(3:end))];
param_file = [param_file, sprintf('Q = %.4e\n', 3e3*npart_save/npart)];
param_file = [param_file, sprintf('npart = %d\n', npart_save)];
param_file = [param_file, sprintf('T = %.6e\n', T)];
param_file = [param_file, sprintf('npt = %d\n', npt)];
param_file = [param_file, sprintf('save_path = %s%s\n', home, datadir(3:end))];
param_file = [param_file, sprintf('save_root = qp_bet\n')];

f = fopen([datadir, 'qp_bet_param.txt'], 'w+');
fwrite(f, param_file);
fclose(f);



%% Compute betatron radiation

system(['time /opt/local/lib/openmpi/bin/mpirun -np ', num2str(n_process), ' ~/Dropbox/SeB/Codes/bin/', bet_executable, ' ', datadir, 'qp_bet_param.txt']);

data = load([datadir, 'qp_bet_ang.txt']);
x = distance_plasma_lanex*unique(data(:,2));
y = distance_plasma_lanex*unique(data(:,3));
dW = reshape(data(:,4), length(x), length(y));
fig = figure(2);
set(fig, 'color', 'w');
set(fig, 'position', [263, 164, 800, 700]);
pcolor(x, y, dW);
colormap(cmap), colorbar(), shading flat;
daspect([1 1 1]);
xlabel('x (mm)'), ylabel('y (mm)');



%% Read sorted beam

if sim_number==0 % for qpic.e.0125_success
    npt = 149;
    npart_save = 125657;
end
if sim_number==19 % for Sim_19
    npt = 1486;
    npart_save = 199649;
end
if sim_number==20 % for Sim_20
    npt = 68;
    npart_save = 262126;
end
if sim_number==21 % for Sim_21
    npt = 149;
    npart_save = 130510;
end
if sim_number==24 % for Sim_24
    npt = 149;
    npart_save = 126777;
end
if sim_number==25 % for Sim_25
    npt = 186;
    npart_save = 129174;
end
if sim_number==26 % for Sim_26
    npt = 186;
    npart_save = 127801;
end
if sim_number==27 % for Sim_27
    npt = 186;
    npart_save = 124192;
end
if sim_number==99 % for Sim_99
    npt = 186;
    npart_save = 131066;
end
if sim_number==102 % for Sim_102
    npt = 425;
    npart_save = 73822;
end

% npt = 66;
% npart_save = 131066;

if npart_save*npt*6 > RAM*1e9/8
    npart_per_step = floor( RAM*1e9/8/(6*npt) );
    n_saving_step = floor(npart_save/npart_per_step);
    reste = npart_save - n_saving_step * npart_per_step;
    SUB_BEAM = zeros(0, npt, 6);
    for i=1:n_saving_step
        istart = 1 + (i-1)*npart_per_step;
        istop = i*npart_per_step;
        BEAM_SORTED = read_qp_traj(datadir, istart, istop, npt);
        tmp = defining_subbeam(BEAM_SORTED);
        disp(size(tmp, 1));
        SUB_BEAM = cat(1, SUB_BEAM, tmp);
    end
    BEAM_SORTED = read_qp_traj(datadir, npart_save-reste+1, npart_save, npt);
    tmp = defining_subbeam(BEAM_SORTED);
    disp(size(tmp, 1));
    SUB_BEAM = cat(1, SUB_BEAM, tmp);
else
    BEAM_SORTED = read_qp_traj(datadir, 1, npart_save, npt);
    SUB_BEAM = defining_subbeam(BEAM_SORTED);
end



%% Get Betatron Profile from SUB_BEAM

datadir2 = [datadir, 'betarad/tmp/'];
SUB_BEAM = defining_subbeam(BEAM_SORTED);

[x, y, dW] = get_betarad(SUB_BEAM, datadir2, home, bet_executable, n_process, npt, T, npart);



%% Get Beam and Wakefield pictures

do_save = 1;
x_range = [-75, 75];

% SI units
scale_E = SI_em*SI_c*omega_p/SI_e;
scale_B = SI_em*SI_c*omega_p/SI_e/SI_c;
Box_X = my_get_quickpic_param(myfile_rpinput, 'Box_X');
Box_Z = my_get_quickpic_param(myfile_rpinput, 'Box_Z');
n_x = 2^(my_get_quickpic_param(myfile_rpinput, 'INDX'));
n_z = 2^(my_get_quickpic_param(myfile_rpinput, 'INDZ'));
% V/m to GV/m
scale_E = scale_E / 1e9;


if do_save
    path_QEB_QEP = [datadir 'movies/QEB-QEP/'];
    mkdir(path_QEB_QEP);
    vidObj_QEB_QEP = VideoWriter([path_QEB_QEP 'QEB-QEP.avi']);
    vidObj_QEB_QEP.FrameRate = 8;
    open(vidObj_QEB_QEP);
    path_EZ_FPERP = [datadir 'movies/EZ-FPERP/'];
    mkdir(path_EZ_FPERP);
    vidObj_EZ_FPERP = VideoWriter([path_EZ_FPERP 'EZ-FPERP.avi']);
    vidObj_EZ_FPERP.FrameRate = 8;
    open(vidObj_EZ_FPERP);
    if potential_on
        path_PHI_PSI = [datadir 'movies/PHI-PSI/'];
        mkdir(path_PHI_PSI);
        vidObj_PHI_PSI = VideoWriter([path_PHI_PSI 'PHI-PSI.avi']);
        vidObj_PHI_PSI.FrameRate = 8;
        open(vidObj_PHI_PSI);
    end
end

waterfall = zeros(n_z,npt-1);

for i=2:npt
    
    s = (i-1)*SI_c*DFPHA_BEAM*DT/omega_p;
    
    n_3D_timestep_str = sprintf('%.4d', n_3D_timestep_start+(i-1)*DFPHA_BEAM);
    if mod(i,10)==0
        disp(n_3D_timestep_str);
    end
    myfile = [datadir 'QEB-XZ/QEB-XZ_' n_3D_timestep_str '.h5'];
    QEB = -double(my_read_hdf(myfile))';
    myfile = [datadir 'QEP1-XZ/QEP1-XZ_' n_3D_timestep_str '.h5'];
    QEP = -double(my_read_hdf(myfile))';
%     myfile = [datadir 'FEX-XZ/FEX-XZ_' n_3D_timestep_str '.h5'];
%     FEX = scale_E * double(my_read_hdf(myfile))';
%     myfile = [datadir 'FEY-XZ/FEY-XZ_' n_3D_timestep_str '.h5'];
%     FEY = scale_E * double(my_read_hdf(myfile))';
    myfile = [datadir 'FEZ-XZ/FEZ-XZ_' n_3D_timestep_str '.h5'];
    FEZ = scale_E * double(my_read_hdf(myfile))';
%     myfile = [datadir 'FBX-XZ/FBX-XZ_' n_3D_timestep_str '.h5'];
%     FBX = scale_B * double(my_read_hdf(myfile))';
%     myfile = [datadir 'FBY-XZ/FBY-XZ_' n_3D_timestep_str '.h5'];
%     FBY = scale_B * double(my_read_hdf(myfile))';
%     myfile = [datadir 'FBZ-XZ/FBZ-XZ_' n_3D_timestep_str '.h5'];
%     FBZ = scale_B * double(my_read_hdf(myfile))';
%     if potential_on
%         myfile = [datadir 'PHI-XZ/PHI-XZ_' n_3D_timestep_str '.h5'];
%         PHI = double(my_read_hdf(myfile))';
%         myfile = [datadir 'PSI-XZ/PSI-XZ_' n_3D_timestep_str '.h5'];
%         PSI = double(my_read_hdf(myfile))';
%     end
%     
    waterfall(:,i-1) = FEZ(n_x/2,:);

    F_perp = -(FEX-1e-9*SI_c*FBY);

    Z = - (1:n_z) * Box_Z/n_z;
    X = (1:n_x) * Box_X/n_x - Box_X/2;
    [ZZ, XX] = meshgrid(Z,X);
    
    fig = figure(1);
    set(fig, 'position', [50, 164, 600, 822]);
    set(fig, 'color', 'w');
    
    
%     if i==2
        subplot(211);
        colormap(wbgyr);
    %     pcolor(ZZ,XX,QEB), shading flat, cb = colorbar();
    %     c_min = min(min((QEB(abs(XX)<x_range(2)))));
    %     c_max = max(max((QEB(abs(XX)<x_range(2)))));
    %     ylim(x_range), caxis([log(c_min), log(c_max)]);
    %     xlabel('z (um)'), ylabel('x (um)'), title(['Beam density (units of n0) at s = ' num2str(1e2*s, '%.2f') ' cm']);
        pcolor(ZZ,XX,log(QEB)), shading flat, cb = colorbar();
%         ax_qeb = get(gca, 'Children');
        ylim(x_range);
        xlabel('z (um)'), ylabel('x (um)'), title(['Beam density (log scale) at s = ' num2str(1e2*s, '%.2f') ' cm']);
        
        subplot(212);
        colormap(wbgyr);
        pcolor(ZZ,XX,QEP), shading flat, cb = colorbar();
%         ax_qep = get(gca, 'Children');
        c_min = min(min((QEP(abs(XX)<x_range(2)))));
        c_max = max(max((QEP(abs(XX)<x_range(2)))));
        caxis([c_min, c_max]);
        ylim(x_range), 
        xlabel('z (um)'), ylabel('x (um)'), title(['Plasma density (units of n0) at s = ' num2str(1e2*s, '%.2f') ' cm']);
%     else
%         set(ax_qeb, 'CData', log(QEB)); 
% %         title(['Beam density (log scale) at s = ' num2str(1e2*s, '%.2f') ' cm']);;
%         set(ax_qep, 'CData', QEP); 
% %         title(['Plasma density (units of n0) at s = ' num2str(1e2*s, '%.2f') ' cm']);
%         c_min = min(min((QEP(abs(XX)<x_range(2)))));
%         c_max = max(max((QEP(abs(XX)<x_range(2)))));
%         caxis([c_min, c_max]);
%     end
    
    if do_save
        writeVideo(vidObj_QEB_QEP, getframe(fig));
        filename = [path_QEB_QEP 'frame_' num2str(i, '%.5d') '.png'];
        saveas(fig, filename, 'png');
    else
        pause(0.001);
    end
%     
%     fig = figure(2);
%     colormap(bwr);
%     set(fig, 'position', [650, 164, 600, 822]);
%     set(fig, 'color', 'w');
% 
%     subplot(211);
%     colormap(bwr);
%     pcolor(ZZ,XX,-FEZ), shading flat, cb = colorbar();
%     c_max = max(max(abs((FEZ(abs(XX)<x_range(2))))));
%     ylim(x_range), caxis([-c_max, c_max]);
%     xlabel('z (um)'), ylabel('x (um)'), title(['F_{long} (GeV/m) at s = ' num2str(1e2*s, '%.2f') ' cm']);
%     
%     subplot(212);
%     colormap(bwr);
%     pcolor(ZZ,XX,F_perp), shading flat, cb = colorbar();
%     c_max = max(max(abs((F_perp(abs(XX)<x_range(2))))));
%     ylim(x_range), caxis([-c_max, c_max]);
%     xlabel('z (um)'), ylabel('x (um)'), title(['F_{perp} (GeV/m) at s = ' num2str(1e2*s, '%.2f') ' cm']);
%    
%     if do_save
%         writeVideo(vidObj_EZ_FPERP, getframe(fig));
%         filename = [path_EZ_FPERP 'frame_' num2str(i, '%.5d') '.png'];
%         saveas(fig, filename, 'png');
%     else
%         pause(0.001);
%     end

%     if potential_on
%         fig = figure(3);
%         colormap(bwr);
%         set(fig, 'position', [1250, 164, 600, 822]);
%         set(fig, 'color', 'w');
%         subplot(211);
%         pcolor(ZZ,XX,PHI), shading flat, colorbar();
%         c_max = max(max(abs((PHI(abs(XX)<x_range(2))))));
%         ylim(x_range), caxis([-c_max, c_max]);
%         xlabel('z (um)'), ylabel('x (um)'), title(['Phi (normalized) at s = ' num2str(1e2*s, '%.2f') ' cm']);
%         subplot(212);
%         pcolor(ZZ,XX,PSI), shading flat, colorbar();
%         c_max = max(max(abs((PSI(abs(XX)<x_range(2))))));
%         ylim(x_range), caxis([-c_max, c_max]);
%         xlabel('z (um)'), ylabel('x (um)'), title(['Psi (normalized) at s = ' num2str(1e2*s, '%.2f') ' cm']);
%         if do_save
%             writeVideo(vidObj_PHI_PSI, getframe(fig));
%             filename = [path_PHI_PSI 'frame_' num2str(i, '%.5d') '.png'];
%             saveas(fig, filename, 'png');
%         else
%             pause(0.001);
%         end
%     end
    
end

if do_save
    close(vidObj_QEB_QEP);
    if potential_on
        close(vidObj_PHI_PSI);
    end
    close(vidObj_EZ_FPERP);
end




%% Waterfall plot of Ez over s

XI = - (1:n_z) * Box_Z/n_z;
S = (1:npt-1) * SI_c*DFPHA_BEAM*DT/omega_p;
[SS, XIXI] = meshgrid(S,XI);

figure(10)
subplot(121)
pcolor(SS, XIXI, waterfall_70cm)
shading interp
caxis([-40, 40]);
cb = colorbar();
colormap(bwr);
cblabel('Ez (GV/m)', 'fontsize', 16)
xlabel('s (m)', 'fontsize', 16)
ylabel('xi (um)', 'fontsize', 16)
title('Beta function of 70 cm', 'fontsize', 16)
axis ij
set(gca, 'fontsize', 16)
hold on;
[C,h] = contour(SS, XIXI, waterfall_70cm, 'k');
subplot(122)
xlabel('s (m)')
ylabel('xi (um)')
pcolor(SS, XIXI, waterfall_5mm)
subplot(122)
shading interp
xlabel('s (m)', 'fontsize', 16)
ylabel('xi (um)', 'fontsize', 16)
title('Beta function of 4.7 mm', 'fontsize', 16)
cb = colorbar();
colormap(bwr);
cblabel('Ez (GV/m)', 'fontsize', 16)
caxis([-40, 40]);
axis ij
set(gca, 'fontsize', 16)
hold on;
[C,h] = contour(SS, XIXI, waterfall_5mm, 'k');


