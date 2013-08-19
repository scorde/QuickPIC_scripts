function [x, y, dW] = get_betarad( SUB_BEAM, datadir, home, bet_executable, n_process, npt, T, Q )

distance_plasma_lanex = 23.27; % in m

mkdir(datadir);

fid = fopen([datadir, 'tmp_traj'], 'wb');
fwrite(fid, permute(SUB_BEAM, [3 2 1]), 'float');
fclose(fid);

npart_save = size(SUB_BEAM, 1);

param_header = sprintf('input_type = 1\nsync_limit = 1\nnpt_rec = 1\nnps = 0\nnps_spacing = 0\ntheta_X_max = 2.5\ntheta_Y_max = 2.5\nangular_resolution = 0.025\nsave_traj = 0\nsave_txt_d2W = 0\n');
param_file = [param_header, sprintf('traj_file = %s%stmp_traj\n', home, datadir(3:end))];
param_file = [param_file, sprintf('Q = %.4e\n', Q)];
param_file = [param_file, sprintf('npart = %d\n', npart_save)];
param_file = [param_file, sprintf('T = %.6e\n', T)];
param_file = [param_file, sprintf('npt = %d\n', npt)];
param_file = [param_file, sprintf('save_path = %s%s\n', home, datadir(3:end))];
param_file = [param_file, sprintf('save_root = tmp_bet\n')];

f = fopen([datadir, 'tmp_bet_param.txt'], 'w+');
fwrite(f, param_file);
fclose(f);

system(['time /opt/local/lib/openmpi/bin/mpirun -np ', num2str(n_process), ' ~/Dropbox/SeB/Codes/bin/', bet_executable, ' ', datadir, 'tmp_bet_param.txt']);

data = load([datadir, 'tmp_bet_ang.txt']);
x = distance_plasma_lanex*unique(data(:,2));
y = distance_plasma_lanex*unique(data(:,3));
dW = reshape(data(:,4), length(x), length(y));


fig = figure(1);
set(fig, 'color', 'w');
pcolor(x, y, dW);
cmap = custom_cmap();
colormap(cmap.wbgyr), colorbar(), shading flat;
daspect([1 1 1]);
xlabel('x (mm)'), ylabel('y (mm)');


system(['rm -R ', datadir]);

end

