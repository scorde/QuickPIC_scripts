%
% Matlab script to analyze QuickPIC output files
% including getting electron trajectories and computing betatron radiation
% 
% S. Corde, Jan 26, 2013
%
% Last update: S. Corde, July 14, 2013
%
% Extraction of electron trajectories and betatron calculation only for QuickPIC binary qpic.e.0125
%



%% Defining input parameters

addpath('~/Dropbox/SeB/Codes/sources/QUICKPICSIM/');

sim_number = 228;
datadir = ['~/QuickPIC_sim/qpic_' num2str(sim_number) '/'];

% Define a memory size to use for storing trajectories, in GB
RAM = 2.3;

% Parameters needed for betatron computation
par_bet = struct();
par_bet.home = '/Users/scorde/';
par_bet.bet_executable = 'bet';
par_bet.n_process = 8;
par_bet.do_spec = 1; % 1 to compute betatron spectrum, 0 otherwise
par_bet.energy_range = [1e4 1e10]; % energy range in eV for betatron spectrum calculation
par_bet.nps = 200; % number of frequency steps for betatron spectrum calculation
par_bet.lanex_response_file = [par_bet.home 'Dropbox/SeB/PostDoc/Projects/BetaRad/scintillator_sim/lanex_sextufilter/lanex_regular_response_W1mm_2x2mm2.txt'];
par_bet.distance_plasma_lanex = 22.2; % in m

% Time step range
npt_start = -1;   % -1: start from start
npt_start = 550;
npt_end = -1;     % -1: go all the way to end
npt_end = 2850;



my_SI_params;
myfile_rpinput = [datadir 'rpinput'];
n0 = my_get_quickpic_param(myfile_rpinput, 'Plasma_Density');
TEND = my_get_quickpic_param(myfile_rpinput, 'TEND');
DT = my_get_quickpic_param(myfile_rpinput, 'DT');
DFPHA_BEAM = my_get_quickpic_param(myfile_rpinput, 'DFPHA_BEAM');
Num_Particle = my_get_quickpic_param(myfile_rpinput, 'Num_Particle');

omega_p = 5.64e4*sqrt(n0);  % Plasma frequency in s-1

if npt_start == -1; npt_start = 0; end;
if npt_end == -1;  npt_end = floor(TEND/DT); end;

npt = floor((npt_end-npt_start)/DFPHA_BEAM) + 1;
T = DFPHA_BEAM*DT*(npt-1)/omega_p;
npart = size(my_get_quickpic_phasespace(datadir, '01', sprintf('%.4d', npt_start+(npt-1)*DFPHA_BEAM)), 1);
npart_tot = size(my_get_quickpic_phasespace(datadir, '01', sprintf('%.4d', 0)), 1);

disp(['Time step range is from ' num2str(npt_start, '%.4d') ' to ' num2str(npt_start+(npt-1)*DFPHA_BEAM, '%.4d')]);



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
        npart_save = npart_save + save_qp_traj(istart, istop, npt_start, npt, datadir, fid);
    end
    npart_save = npart_save + save_qp_traj(npart-reste+1, npart, npt_start, npt, datadir, fid);
    fclose(fid);
else
    fid = fopen([datadir 'qp_traj'], 'wb');
    npart_save = save_qp_traj(1, npart, npt_start, npt, datadir, fid);
    fclose(fid);
end

Q = Num_Particle * SI_e * 1e12 *npart_save/npart_tot;

disp('Total number of macro-particles saved:');
disp(npart_save);
disp('Number of time steps saved:');
disp(npt);
disp('Total time T to input in the betatron code');
disp(T);
disp('Beam charge Q to input in the betatron code');
disp(Q);

save([datadir 'qp_traj_param.mat'], 'npart', 'npart_save', 'npart_tot', 'npt_start', 'npt', 'T', 'Q');



%% Compute betatron radiation

par_bet.datadir = datadir;
par_bet.datadir_bet = [datadir 'betarad_' num2str(npt_start, '%.4d') '_' num2str(npt_end, '%.4d') '/'];
par_bet.rm_dir = 0;
get_betarad(par_bet);



%% Read sorted beam and select a SUB_BEAM

par_subbeam = struct();
par_subbeam.energy = 10;
par_subbeam.dE = 20;
par_subbeam.slice = 0;
par_subbeam.ds = 1000;

load([datadir 'qp_traj_param.mat']);

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
    tmp = defining_subbeam(BEAM_SORTED, par_subbeam);
    disp(size(tmp, 1));
    SUB_BEAM = cat(1, SUB_BEAM, tmp);
else
    BEAM_SORTED = read_qp_traj(datadir, 1, npart_save, npt);
    SUB_BEAM = defining_subbeam(BEAM_SORTED, par_subbeam);
end



%% Select a SUB_BEAM

par_subbeam = struct();
par_subbeam.energy = 10;
par_subbeam.dE = 20;
par_subbeam.slice = 0;
par_subbeam.ds = 1000;
SUB_BEAM = defining_subbeam(BEAM_SORTED, par_subbeam);

    

%% Compute betatron radiation from SUB_BEAM

par_bet.datadir = datadir;
par_bet.datadir_bet = [datadir 'betarad_no-acc_' num2str(npt_start, '%.4d') '_' num2str(npt_end, '%.4d') '/'];
par_bet.Q = Num_Particle * SI_e * 1e12 * size(SUB_BEAM,1)/npart_tot;
par_bet.T = T;
par_bet.rm_dir = 0;
betarad = get_betarad(par_bet, SUB_BEAM);



%% Get Beam and Wakefield pictures

par = custom_cmap();

par.do_save = 1;
do_log_QEB = 1;
do_log_QEP = 0;
par.x_range = [-150, 150];
par.fontsize = 14;


if my_get_quickpic_param(myfile_rpinput, 'EX0') == 0; par.slice = 'XZ'; end;
if my_get_quickpic_param(myfile_rpinput, 'EY0') == 0; par.slice = 'YZ'; end;

figure(1);
set(1, 'position', [50, 164, 600, 822]);
set(1, 'PaperPosition', [0., 0., 6, 8]);
set(1, 'color', 'w');
clf();
figure(2);
set(2, 'position', [650, 164, 600, 822]);
set(2, 'PaperPosition', [0., 0., 6, 8]);
set(2, 'color', 'w');

% SI units
scale_E = SI_em*SI_c*omega_p/SI_e;
scale_B = SI_em*SI_c*omega_p/SI_e/SI_c;
Box_X = my_get_quickpic_param(myfile_rpinput, ['Box_' par.slice(1)]);
Box_Z = my_get_quickpic_param(myfile_rpinput, 'Box_Z');
n_x = 2^(my_get_quickpic_param(myfile_rpinput, ['IND'  par.slice(1)]));
n_z = 2^(my_get_quickpic_param(myfile_rpinput, 'INDZ'));
% V/m to GV/m
scale_E = scale_E / 1e9;

[par.ZZ, par.XX] = meshgrid(- (1:n_z) * Box_Z/n_z, (1:n_x) * Box_X/n_x - Box_X/2);

waterfall = zeros(n_z,npt-1);

path_QEB_QEP = [datadir 'movies/QEB-QEP/'];
mkdir(path_QEB_QEP);
path_EZ_FPERP = [datadir 'movies/EZ-FPERP/'];
mkdir(path_EZ_FPERP);


for i=2:npt
    
    s = npt_start*SI_c*DT/omega_p + (i-1)*SI_c*DFPHA_BEAM*DT/omega_p;
    s_str = ['s = ' num2str(1e2*s, '%.2f') ' cm  '];
    
    npt_str = sprintf('%.4d', npt_start+(i-1)*DFPHA_BEAM);
    if mod(i,10)==0; disp(npt_str); end;
    
    myfile = [datadir 'QEB-' par.slice '/QEB-' par.slice '_' npt_str '.h5'];
    QEB = -double(my_read_hdf(myfile))';
    myfile = [datadir 'QEP1-' par.slice '/QEP1-' par.slice '_' npt_str '.h5'];
    QEP = -double(my_read_hdf(myfile))';
    myfile = [datadir 'FEX-' par.slice '/FEX-' par.slice '_' npt_str '.h5'];
    FEX = scale_E * double(my_read_hdf(myfile))';
    myfile = [datadir 'FEY-' par.slice '/FEY-' par.slice '_' npt_str '.h5'];
    FEY = scale_E * double(my_read_hdf(myfile))';
    myfile = [datadir 'FEZ-' par.slice '/FEZ-' par.slice '_' npt_str '.h5'];
    FEZ = scale_E * double(my_read_hdf(myfile))';
    myfile = [datadir 'FBX-' par.slice '/FBX-' par.slice '_' npt_str '.h5'];
    FBX = scale_B * double(my_read_hdf(myfile))';
    myfile = [datadir 'FBY-' par.slice '/FBY-' par.slice '_' npt_str '.h5'];
    FBY = scale_B * double(my_read_hdf(myfile))';
    myfile = [datadir 'FBZ-' par.slice '/FBZ-' par.slice '_' npt_str '.h5'];
    FBZ = scale_B * double(my_read_hdf(myfile))';

    waterfall(:,i-1) = FEZ(n_x/2,:);

    if strcmp(par.slice, 'XZ'); F_perp = -(FEX-1e-9*SI_c*FBY); end;
    if strcmp(par.slice, 'YZ'); F_perp = -(FEY+1e-9*SI_c*FBX); end;
        
    % Plot QEB and QEP
    par.path = path_QEB_QEP;
    par.title = {{['  Ez (units of n0) at ' s_str], ['  Beam density (log scale) at ' s_str]},...
        {['  Plasma density (units of n0) at ' s_str], ['  Plasma density (log scale) at ' s_str]}};
    par.fig = 1;
    if i==2
        ax_struct_QEB_QEP = dual_plot_ini(par, QEB, QEP, do_log_QEB, do_log_QEP, i, 'wbgyr');
    else
        dual_plot_set(ax_struct_QEB_QEP, par, QEB, QEP, do_log_QEB, do_log_QEP, i, 'wbgyr');
    end
    
    % Plot F_long and F_perp
    par.path = path_EZ_FPERP;
    par.title = {{['  F_{long} (GeV/m) at ' s_str], ['  ' s_str]},...
        {['  F_{perp} (GeV/m) at ' s_str], ['  ' s_str]}};
    par.fig = 2;
    if i==2
        ax_struct_EZ_FPERP = dual_plot_ini(par, -FEZ, F_perp, 0, 0, i, 'bwr');
    else
        dual_plot_set(ax_struct_EZ_FPERP, par, -FEZ, F_perp, 0, 0, i, 'bwr');
    end

end

if par.do_save
    system(['/usr/local/bin/ffmpeg -i ' path_QEB_QEP 'frame_%05d.png ' datadir 'movies/qpic_' num2str(sim_number) '_QEB_QEP.mpeg']);
    system(['/usr/local/bin/ffmpeg -i ' path_EZ_FPERP 'frame_%05d.png ' datadir 'movies/qpic_' num2str(sim_number) '_EZ_FPERP.mpeg']);
end



%% Waterfall plot of Ez over s

XI = - (1:n_z) * Box_Z/n_z;
S = (1:npt-1) * SI_c*DFPHA_BEAM*DT/omega_p;
[SS, XIXI] = meshgrid(S,XI);

fig = figure(4);
set(fig, 'color', 'w');
set(fig, 'position', [263, 164, 800, 700]);
pcolor(SS, XIXI, waterfall);
shading interp;
caxis([-50, 50]);
cb = colorbar();
colormap(par.bwr);
% cblabel('Ez (GV/m)', 'fontsize', 16);
xlabel('s (m)', 'fontsize', 16);
ylabel('xi (um)', 'fontsize', 16);
title('Waterfall plot of on-axis E_z', 'fontsize', 16);
axis ij;
set(gca, 'fontsize', 16);
hold on;
contour(SS, XIXI, waterfall, 'k');



%% Get longitudinal phase space pictures

do_save = 1;

xi_axis = linspace(-70, 70);
E = linspace(0, 60);

fig = figure(5);
set(fig, 'color', 'w');
set(fig, 'position', [263, 164, 800, 700]);
set(fig, 'PaperPosition', [0., 0., 8, 8]);
set(gca, 'fontsize', 16);
mkdir([datadir 'movies/longitudinal/']);

for i=1:size(BEAM_SORTED, 2)
	s = npt_start*SI_c*DT/omega_p + (i-1)*SI_c*DFPHA_BEAM*DT/omega_p;
%     xi = SUB_BEAM(:,i,3) - 1e6*(i-1)*SI_c*DFPHA_BEAM*DT/omega_p;
%     histmat = hist2(xi, SUB_BEAM(:,i,6), xi_axis, E);
    histmat = hist2(BEAM_SORTED(:,1,3), BEAM_SORTED(:,i,6), xi_axis, E);
    pcolor(xi_axis(2:end), E(2:end), log10(histmat(2:end, 2:end)));
    cmap = custom_cmap();
    colormap(cmap.wbgyr), shading flat;
    xlabel('z (um)  '), ylabel('E (GeV)  '), title(['Longitudinal phase space at s = ' num2str(1e2*s, '%.2f') ' cm  ']);
    if do_save
        filename = [datadir 'movies/longitudinal/frame_' num2str(i, '%.5d') '.png'];
        saveas(fig, filename, 'png');
    else
        pause(0.001);
    end
end

if do_save
    system(['/usr/local/bin/ffmpeg -i ' datadir 'movies/longitudinal/frame_%05d.png ' datadir 'movies/qpic_' num2str(sim_number) '_longitudinal.mpeg']);
end








