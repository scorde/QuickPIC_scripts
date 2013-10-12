function output = get_betarad( par_bet, SUB_BEAM )


fontsize = 18;
cmap = custom_cmap();

mkdir(par_bet.datadir_bet);

if nargin>1
    fid = fopen([par_bet.datadir_bet, 'qp_traj'], 'wb');
    fwrite(fid, permute(SUB_BEAM, [3 2 1]), 'float');
    fclose(fid);
    par_bet.npart_save = size(SUB_BEAM, 1);
    par_bet.npt = size(SUB_BEAM, 2);
else
    load([par_bet.datadir 'qp_traj_param.mat']);
    par_bet.Q = Q;
    par_bet.npart_save = npart_save;
    par_bet.T = T;
    par_bet.npt = npt;
end


param_header = sprintf(['input_type = 1\nsync_limit = 1\nnpt_rec = 1'...
    '\nnps_spacing = 0\ntheta_X_max = 2.5\ntheta_Y_max = 2.5'...
    '\nangular_resolution = 0.025\nsave_traj = 0\nsave_txt_d2W = 0\n']);
if nargin>1
    param_file = [param_header, sprintf('traj_file = %s%sqp_traj\n', par_bet.home, par_bet.datadir_bet(3:end))];
else
    param_file = [param_header, sprintf('traj_file = %s%sqp_traj\n', par_bet.home, par_bet.datadir(3:end))];
end
param_file = [param_file, sprintf('Q = %.4e\n', par_bet.Q)];
param_file = [param_file, sprintf('npart = %d\n', par_bet.npart_save)];
param_file = [param_file, sprintf('T = %.6e\n', par_bet.T)];
param_file = [param_file, sprintf('npt = %d\n', par_bet.npt)];
param_file = [param_file, sprintf('omega_min = %.4e\n', par_bet.energy_range(1))];
param_file = [param_file, sprintf('omega_max = %.4e\n', par_bet.energy_range(2))];
if par_bet.do_spec
    param_file = [param_file, sprintf('nps = %d\n', par_bet.nps)];
else
    param_file = [param_file, sprintf('nps = %d\n', 0)];
end
param_file = [param_file, sprintf('save_path = %s%s\n', par_bet.home, par_bet.datadir_bet(3:end))];
param_file = [param_file, sprintf('save_root = qp_bet\n')];

f = fopen([par_bet.datadir_bet, 'qp_bet_param.txt'], 'w+');
fwrite(f, param_file);
fclose(f);

system(['time /opt/local/lib/openmpi/bin/mpirun -np ', num2str(par_bet.n_process), ' ~/Dropbox/SeB/Codes/bin/', par_bet.bet_executable, ' ', par_bet.datadir_bet, 'qp_bet_param.txt']);

data = load([par_bet.datadir_bet, 'qp_bet_ang.txt']);
x = par_bet.distance_plasma_lanex*unique(data(:,2));
y = par_bet.distance_plasma_lanex*unique(data(:,3));
dW = reshape(data(:,4), length(y), length(x));
fig = figure(3);
set(fig, 'color', 'w');
set(fig, 'position', [263, 164, 800, 700]);
set(fig, 'PaperPosition', [0., 0., 8, 8]);
set(gca, 'fontsize', fontsize);
pcolor(x, y, dW);
colormap(cmap.wbgyr), colorbar('fontsize', fontsize), shading flat;
daspect([1 1 1]);
xlabel('x (mm)'), ylabel('y (mm)');
caxis([0 max(dW(:))]);
saveas(fig, [par_bet.datadir_bet 'qp_bet_ang'], 'png');

if par_bet.do_spec
    E = logspace(log10(par_bet.energy_range(1)), log10(par_bet.energy_range(2)), par_bet.nps);
    tmp = load(par_bet.lanex_response_file);
    lanex_response = interp1(tmp(:,1), tmp(:,4), 1e-6*E)';
    fid = fopen([par_bet.datadir_bet 'qp_bet_d2W'], 'rb');
    d2W = fread(fid, par_bet.nps*(201)^2, 'double');
    fclose(fid);
    d2W = reshape(d2W, 201, 201, par_bet.nps);
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
    saveas(fig, [par_bet.datadir_bet 'qp_bet_dS'], 'png');
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
    saveas(fig, [par_bet.datadir_bet 'qp_bet_dS_sym'], 'png');
    dS_filt = filter2(ones(5,5)/25, dS);
    qp_bet = struct();
    qp_bet.gamma_max_1 = max(dS_filt(:));
    tmp_1 = dS_filt>qp_bet.gamma_max_1/2.;
    gamma_div_1 = 2*sqrt( sum(tmp_1(:))*(x(2)-x(1))*(y(2)-y(1))/pi );
    qp_bet.gamma_yield_1 = qp_bet.gamma_max_1*gamma_div_1^2;
    qp_bet.gamma_div_1 = gamma_div_1/par_bet.distance_plasma_lanex;
    qp_bet.gamma_max_2 = max(dS_sym(:));
    tmp_2 = dS_sym>qp_bet.gamma_max_2/2.;
    gamma_div_2 = 2*sqrt( sum(tmp_2(:))*(x(2)-x(1))*(y(2)-y(1))/pi );
    qp_bet.gamma_yield_2 = qp_bet.gamma_max_2*gamma_div_2^2;
    qp_bet.gamma_div_2 = gamma_div_2/par_bet.distance_plasma_lanex;
    disp(qp_bet);
    save([par_bet.datadir_bet 'qp_bet.mat'], 'qp_bet');
end



if nargin>1; system(['rm ', par_bet.datadir_bet, 'qp_traj']); end;
if par_bet.rm_dir; system(['rm -R ', par_bet.datadir_bet]); end;


output.x = x;
output.y = y;
output.E = E;
output.dW = dW;
output.d2W = d2W;
output.qp_bet = qp_bet;


end

