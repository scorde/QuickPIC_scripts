function SUB_BEAM = defining_subbeam( BEAM_SORTED, par_in )

par.energy = 25;
par.dE = 2.;

par.slice = 0.;
par.ds = 1000;

% cond = abs(BEAM_SORTED(:,end,6))>20+dE/2.;
% cond = abs(BEAM_SORTED(:,end,6))<20-dE/2.;
% cond = abs(20-BEAM_SORTED(:,end,6))>dE/2.;

if nargin > 1; par = par_in; end;

cond = abs(par.energy-BEAM_SORTED(:,end,6))<par.dE/2. & abs(par.slice-BEAM_SORTED(:,1,3))<par.ds/2.;

SUB_BEAM = BEAM_SORTED(cond, :, :);

end

