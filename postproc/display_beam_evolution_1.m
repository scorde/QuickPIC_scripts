


%% Display evolution of beam parameters

do_save = 0;
save_name = '';

energy = linspace(22.5, 22.5, 1);
dE = 500;

slice = linspace(0, 0, 1);
ds = 3;


Z = (1:size(BEAM_SORTED,2))*SI_c*DFPHA_BEAM*DT/omega_p;

for i=energy
    for j=slice
%         cond = BEAM_SORTED(:,end,6)>20+0.5;
        cond = abs(i-BEAM_SORTED(:,end,6))<dE/2. & abs(j-BEAM_SORTED(:,1,3))<ds/2.;
        
        SUB_BEAM = BEAM_SORTED(cond,:,:);
        beam_param = get_beam_param( SUB_BEAM, size(BEAM_SORTED,1) );
        disp([beam_param.Fraction]);
        
        fig = figure(1);
        set(fig, 'position', [263, 164, 1395, 822]);    
        subplot(231)
        plot(Z, beam_param.sig_x, 'b'), hold on
        plot(Z, beam_param.sig_y, 'r');
        plot([Z(1), Z(end)], [40, 40], 'g-'), hold off
        set(gca, 'fontsize', 20);
        legend('x','y'), xlabel('z (m)'), ylabel('Beam size (um)');
        xlim([0, 0.3]);%, ylim([0, 50]);
        subplot(232)
        plot(Z, beam_param.sig_xp, 'b'), hold on
        plot(Z, beam_param.sig_yp, 'r');
        plot([Z(1), Z(end)], [127.75, 127.75], 'g-');
        plot([Z(1), Z(end)], [383.25, 383.25], 'g-'), hold off
        set(gca, 'fontsize', 20);
        legend('x','y'), xlabel('z (m)'), ylabel('Angular spread (urad)');
        xlim([0, 0.3]);
        subplot(233)
        plot(Z, 1e-6*beam_param.eps_nx, 'b'), hold on
        plot(Z, 1e-6*beam_param.eps_ny, 'r');
        plot([Z(1), Z(end)], [200, 200], 'g-');
        plot([Z(1), Z(end)], [600, 600], 'g-'), hold off
        set(gca, 'fontsize', 20);
        legend('x','y'), xlabel('z (m)'), ylabel('Normalized emittance (um)');
        xlim([0, 0.3]);
        subplot(234)
        plot(Z, beam_param.beta_x, 'b'), hold on
        plot(Z, beam_param.beta_y, 'r');
        plot([Z(1), Z(end)], [0.1044, 0.1044], 'g-');
        plot([Z(1), Z(end)], [0.3131, 0.3131], 'g-'), hold off
        set(gca, 'fontsize', 20);
        legend('x','y'), xlabel('z (m)'), ylabel('Beta function (m)');
        xlim([0, 0.3]);
        subplot(235)
        plot(Z, beam_param.alpha_x, 'b'), hold on
        plot(Z, beam_param.alpha_y, 'r'), hold off
        set(gca, 'fontsize', 20);
        legend('x','y'), xlabel('z (m)'), ylabel('Alpha');
        xlim([0, 0.3]);
        pause;
    end
end






