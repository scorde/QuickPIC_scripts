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
addpath('~/Dropbox/SeB/Codes/sources/E200_scripts/tools/');
addpath('~/Dropbox/SeB/Codes/sources/QuickPIC_scripts/postproc/');

sim_number = 220;
datadir = ['~/QuickPIC_sim/qpic_' num2str(sim_number) '/'];

% Define a memory size to use for storing trajectories, in GB
RAM = 2.5;

% Parameters needed for betatron computation
home = '/Users/scorde/';
bet_executable = 'bet';
n_process = 8;
do_spec = 1; % 1 to compute betatron spectrum, 0 otherwise
energy_range = [1e4 1e10]; % energy range in eV for betatron spectrum calculation
nps = 200; % number of frequency steps for betatron spectrum calculation
lanex_response_file = [home 'Dropbox/SeB/PostDoc/Projects/BetaRad/scintillator_sim/lanex_sextufilter/lanex_regular_response_W1mm_2x2mm2.txt'];
distance_plasma_lanex = 22.2; % in m

% Time step range
npt_start = -1;   % -1: start from start
% npt_start = 1000;
npt_end = -1;     % -1: go all the way to end
% npt_end = 6000;



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

Q = Num_Particle * SI_e * 1e12 *npart_save/npart;

disp('Total number of macro-particles saved:');
disp(npart_save);
disp('Number of time steps saved:');
disp(npt);
disp('Total time T to input in the betatron code');
disp(T);
disp('Beam charge Q to input in the betatron code');
disp(Q);

save([datadir 'qp_traj_param.mat'], 'npart', 'npart_save', 'npt_start', 'npt', 'T', 'Q');



%% Compute betatron radiation

fontsize = 18;

load([datadir 'qp_traj_param.mat']);

param_header = sprintf(['input_type = 1\nsync_limit = 1\nnpt_rec = 1'...
    '\nnps_spacing = 0\ntheta_X_max = 2.5\ntheta_Y_max = 2.5'...
    '\nangular_resolution = 0.025\nsave_traj = 0\nsave_txt_d2W = 0\n']);
param_file = [param_header, sprintf('traj_file = %s%sqp_traj\n', home, datadir(3:end))];
param_file = [param_file, sprintf('Q = %.4e\n', Q)];
param_file = [param_file, sprintf('npart = %d\n', npart_save)];
param_file = [param_file, sprintf('T = %.6e\n', T)];
param_file = [param_file, sprintf('npt = %d\n', npt)];
param_file = [param_file, sprintf('omega_min = %.4e\n', energy_range(1))];
param_file = [param_file, sprintf('omega_max = %.4e\n', energy_range(2))];
if do_spec
    param_file = [param_file, sprintf('nps = %d\n', nps)];
else
    param_file = [param_file, sprintf('nps = %d\n', 0)];
end
param_file = [param_file, sprintf('save_path = %s%sbetarad/\n', home, datadir(3:end))];
param_file = [param_file, sprintf('save_root = qp_bet\n')];

mkdir([datadir 'betarad']);
f = fopen([datadir, 'betarad/qp_bet_param.txt'], 'w+');
fwrite(f, param_file);
fclose(f);

system(['time /opt/local/lib/openmpi/bin/mpirun -np ', num2str(n_process), ' ~/Dropbox/SeB/Codes/bin/', bet_executable, ' ', datadir, 'betarad/qp_bet_param.txt']);

data = load([datadir, 'betarad/qp_bet_ang.txt']);
x = distance_plasma_lanex*unique(data(:,2));
y = distance_plasma_lanex*unique(data(:,3));
dW = reshape(data(:,4), length(y), length(x));
fig = figure(3);
set(fig, 'color', 'w');
set(fig, 'position', [263, 164, 800, 700]);
set(fig, 'PaperPosition', [0., 0., 8, 8]);
set(gca, 'fontsize', fontsize);
pcolor(x, y, dW);
cmap = custom_cmap();
colormap(cmap.wbgyr), colorbar('fontsize', fontsize), shading flat;
daspect([1 1 1]);
xlabel('x (mm)'), ylabel('y (mm)');
caxis([0 max(dW(:))]);
saveas(fig, [datadir 'betarad/qp_bet_ang'], 'png');

if do_spec
    E = logspace(log10(energy_range(1)), log10(energy_range(2)), nps);
    tmp = load(lanex_response_file);
    lanex_response = interp1(tmp(:,1), tmp(:,4), 1e-6*E)';
    fid = fopen([datadir 'betarad/qp_bet_d2W'], 'rb');
    d2W = fread(fid, nps*(201)^2, 'double');
    fclose(fid);
    d2W = reshape(d2W, 201, 201, nps);
    d2W = permute(d2W, [3 1 2]);
    dS = squeeze(trapz(E,d2W.*repmat(lanex_response, [1 201 201]),1));
    dS_sym = zeros(size(dS));
    for i = 0:359; dS_sym = dS_sym + (1/360)*imrotate_sc(dS, i); end;
    fig = figure(6);
    set(fig, 'color', 'w');
    set(fig, 'position', [263, 164, 800, 700]);
    set(fig, 'PaperPosition', [0., 0., 8, 8]);
    set(gca, 'fontsize', fontsize);
    pcolor(x, y, dS);
    colormap(cmap.wbgyr), colorbar('fontsize', fontsize), shading flat;
    daspect([1 1 1]);
    xlabel('x (mm)'), ylabel('y (mm)');
    caxis([0 max(dS(:))]);
    saveas(fig, [datadir 'betarad/qp_bet_dS'], 'png');
    fig = figure(7);
    set(fig, 'color', 'w');
    set(fig, 'position', [263, 164, 800, 700]);
    set(fig, 'PaperPosition', [0., 0., 8, 8]);
    set(gca, 'fontsize', fontsize);
    pcolor(x, y, dS_sym);
    colormap(cmap.wbgyr), colorbar('fontsize', fontsize), shading flat;
    daspect([1 1 1]);
    xlabel('x (mm)'), ylabel('y (mm)');
    caxis([0 max(dS_sym(:))]);
    saveas(fig, [datadir 'betarad/qp_bet_dS_sym'], 'png');
    qp_bet = struct();
    qp_bet.gamma_max = max(dS_sym(:));
    tmp = dS_sym>qp_bet.gamma_max/2.;
    gamma_div = 2*sqrt( sum(tmp(:))*(x(2)-x(1))*(y(2)-y(1))/pi );
    qp_bet.gamma_yield = qp_bet.gamma_max*gamma_div^2;
    qp_bet.gamma_div = gamma_div/distance_plasma_lanex;
    disp(qp_bet);
    save([datadir 'betarad/qp_bet.mat'], 'qp_bet');
end



%% Read sorted beam and select a SUB_BEAM

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
    tmp = defining_subbeam(BEAM_SORTED);
    disp(size(tmp, 1));
    SUB_BEAM = cat(1, SUB_BEAM, tmp);
else
    BEAM_SORTED = read_qp_traj(datadir, 1, npart_save, npt);
    SUB_BEAM = defining_subbeam(BEAM_SORTED);
end



%% Select a SUB_BEAM

par = struct();
par.energy = 30;
par.dE = 2;
par.slice = 0;
par.ds = 1000;
SUB_BEAM = defining_subbeam(BEAM_SORTED, par);

    

%% Get Betatron Profile from SUB_BEAM

[x, y, dW] = get_betarad(SUB_BEAM, [datadir, 'betarad/tmp/'], home, bet_executable,...
    n_process, npt, T, npart);



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








