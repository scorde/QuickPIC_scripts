


%% Choose a subset of the beam

cond = (BEAM_SORTED(:,end,1).^2+BEAM_SORTED(:,end,2).^2)<100 & abs(BEAM_SORTED(:,end,3)-55)<3.;

cond1 = abs(20-BEAM_SORTED(:,end,6))>1;

cond2 = max(abs(BEAM_SORTED(:,:,1)),[],2) < 400;
cond3 = max(abs(BEAM_SORTED(:,:,2)),[],2) < 400;

% cond = cond1 & cond2 & cond3;
cond = cond2 & cond3;

SUB_BEAM = BEAM_SORTED;
% SUB_BEAM = BEAM_SORTED(cond,:,:);



%% Display traj

while 1
    i = randperm(size(SUB_BEAM,1),1);
    plot(SUB_BEAM(i,:,3), SUB_BEAM(i,:,1));
    hold on;
    plot([0, 3e5], [-450, -450], 'r-');
    plot([0, 3e5], [450, 450], 'r-');
%     plot3(SUB_BEAM(i,:,3), SUB_BEAM(i,:,1), SUB_BEAM(i,:,2));
    xlim([-0.1e4 30e4]), ylim([-500 500]), zlim([-500 500]);
%     xlim([-0.1e4 30e4]), ylim([-100 100]), zlim([-100 100]);
    grid on;
%     daspect([3e2 1 1]);
    hold off;
    pause(0.5);
end



