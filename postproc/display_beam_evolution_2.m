


%% Display longitudinal evolution

do_save = 1;

xi_axis = linspace(-150, 150);
E = linspace(0, 40);
fig = figure(5);
set(fig, 'color', 'w');
set(gca, 'fontsize', 20);
if do_save
    mkdir([datadir 'movies/longitudinal/']);
    vidObj = VideoWriter([datadir 'movies/longitudinal/longitudinal.avi']);
    vidObj.FrameRate = 10;
    open(vidObj);
end
for i=1:size(BEAM_SORTED, 2)
	s = (i-1)*SI_c*DFPHA_BEAM*DT/omega_p;
%     xi = SUB_BEAM(:,i,3) - 1e6*(i-1)*SI_c*DFPHA_BEAM*DT/omega_p;
%     histmat = hist2(xi, SUB_BEAM(:,i,6), xi_axis, E);
    histmat = hist2(BEAM_SORTED(:,1,3), BEAM_SORTED(:,i,6), xi_axis, E);
    pcolor(xi_axis(2:end), E(2:end), log10(histmat(2:end, 2:end)));
    colormap(cmap), shading flat;
    xlabel('z (um)'), ylabel('E (GeV)'), title(['Longitudinal phase space at s = ' num2str(1e2*s, '%.2f') ' cm']);
    if do_save
        writeVideo(vidObj,getframe(fig));
        filename = [datadir 'movies/longitudinal/frame_' num2str(i, '%.5d') '.png'];
        saveas(fig, filename, 'png');
    else
        pause(0.001);
    end
end
if do_save
    close(vidObj);
end



%% Display transverse evolution

do_save = 1;
save_name = 's=-30_E=all';

energy = 22.5;
dE = 100.;

slice = -30.;
ds = 2.;


% cond = abs(BEAM_SORTED(:,end,6))>20+dE/2.;
% cond = abs(BEAM_SORTED(:,end,6))<20-dE/2.;
% cond = abs(20-BEAM_SORTED(:,end,6))>dE/2.;
cond = abs(energy-BEAM_SORTED(:,end,6))<dE/2. & abs(slice-BEAM_SORTED(:,1,3))<ds/2.;

SUB_BEAM = BEAM_SORTED(cond,:,:);


xi = linspace(-150, 150);
E = linspace(15, 28);
x = linspace(-50, 50);
y = linspace(-50, 50);
xp = linspace(-4e3, 4e3);
yp = linspace(-4e3, 4e3);
r_axis = linspace(0,120);
rp_axis = linspace(0,4e3);
r_init = sqrt(SUB_BEAM(:,1,1).^2+SUB_BEAM(:,1,2).^2);
Lz_axis = linspace(-5e3,5e3);

fig = figure(2);
set(fig, 'position', [263, 164, 1395, 822]);
% set(fig, 'OuterPosition', [263, 164, 1395, 822]);
set(fig, 'color', 'w');
if do_save
    path = [datadir 'movies/' save_name '/'];
    mkdir(path);
    vidObj = VideoWriter([path save_name '.avi']);
    vidObj.FrameRate = 8;
    open(vidObj);
end
for i=1:size(BEAM_SORTED, 2)
    subplot(241), set(gca, 'fontsize', 20);
    subplot(242), set(gca, 'fontsize', 20);
    subplot(243), set(gca, 'fontsize', 20);
    subplot(244), set(gca, 'fontsize', 20);
    subplot(245), set(gca, 'fontsize', 20);
    subplot(246), set(gca, 'fontsize', 20);
    subplot(247), set(gca, 'fontsize', 20);
    subplot(248), set(gca, 'fontsize', 20);

    r = sqrt(SUB_BEAM(:,i,1).^2+SUB_BEAM(:,i,2).^2);
    rp = sqrt(SUB_BEAM(:,i,4).^2+SUB_BEAM(:,i,5).^2);
    Lz = 1e-3/0.511 * SUB_BEAM(:,i,6)  .* (SUB_BEAM(:,i,1).*SUB_BEAM(:,i,5) - SUB_BEAM(:,i,2).*SUB_BEAM(:,i,4));

    subplot(241);
    histmat = hist2(SUB_BEAM(:,i,1), SUB_BEAM(:,i,2), x, y);
    pcolor(x(2:end), y(2:end), log10(histmat(2:end, 2:end)));
    daspect([1 1 1]);
    colormap(cmap), shading flat;
    xlabel('x (um)'), ylabel('y (um)');
    subplot(242);
    histmat = hist2(SUB_BEAM(:,i,4), SUB_BEAM(:,i,5), xp, yp);
    pcolor(xp(2:end), yp(2:end), log10(histmat(2:end, 2:end)));
    daspect([1 1 1]);
    colormap(cmap), shading flat;
    xlabel('xp (urad)'), ylabel('yp (urad)');
    subplot(243);
    histmat = hist2(SUB_BEAM(:,1,3), SUB_BEAM(:,i,6), xi, E);
    pcolor(xi(2:end), E(2:end), log10(histmat(2:end, 2:end)));
    colormap(cmap), shading flat;
    xlabel('z (um)'), ylabel('E (GeV)');
    subplot(244);
    histmat = hist2(r_init, SUB_BEAM(:,i,6), r_axis, E);
    pcolor(r_axis(2:end), E(2:end),log10(histmat(2:end, 2:end)));
    colormap(cmap), shading flat;
    xlabel('r\_init (um)'), ylabel('E (GeV)');  
    subplot(245);
    histmat = hist2(SUB_BEAM(:,i,1), SUB_BEAM(:,i,4), x, xp);
    pcolor(x(2:end), xp(2:end), log10(histmat(2:end, 2:end)));
    colormap(cmap), shading flat;
    xlabel('x (um)'), ylabel('xp (urad)');
    subplot(246);
    histmat = hist2(SUB_BEAM(:,i,2), SUB_BEAM(:,i,5), y, yp);
    pcolor(y(2:end), yp(2:end), log10(histmat(2:end, 2:end)));
    colormap(cmap), shading flat;
    xlabel('y (um)'), ylabel('yp (urad)');
    subplot(247);
    histmat = hist2(r, rp, r_axis, rp_axis);
    pcolor(r_axis(2:end), rp_axis(2:end), log10(histmat(2:end, 2:end)));
%     plot(r_axis(2:end), sum(histmat(2:end, 2:end),1));
    colormap(cmap), shading flat;
    xlabel('r (um)'), ylabel('rp (urad)');
    subplot(248);
    histmat = hist(Lz, Lz_axis);
    plot(Lz_axis(2:end-1), log10(histmat(2:end-1)));
    xlabel('Lz (um)'), ylabel('dN/dLz (arb. u.)');
    
    if do_save
        writeVideo(vidObj, getframe(fig));
        filename = [path 'frame_' num2str(i, '%.5d') '.png'];
        saveas(fig, filename, 'png');
        print('-dpng', filename);
    else
        pause(0.001);
    end
end
if do_save
    close(vidObj);
end









