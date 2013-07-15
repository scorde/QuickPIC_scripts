function SUB_BEAM = defining_subbeam( BEAM_SORTED )

energy = 17.5;
dE = 100.;

slice = 0.;
ds = 10;

% cond = abs(BEAM_SORTED(:,end,6))>20+dE/2.;
% cond = abs(BEAM_SORTED(:,end,6))<20-dE/2.;
% cond = abs(20-BEAM_SORTED(:,end,6))>dE/2.;
cond = abs(energy-BEAM_SORTED(:,end,6))<dE/2. & abs(slice-BEAM_SORTED(:,1,3))<ds/2.;

SUB_BEAM = BEAM_SORTED(cond, :, :);

end

