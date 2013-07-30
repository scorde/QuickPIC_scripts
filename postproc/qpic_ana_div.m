



par = struct();
par.energy = 30;
par.dE = 1;
par.slice = 0;
par.ds = 1000;

E = [];
div_x = [];
div_y = [];
for i = 1:40
    E(end+1) = i;
    par.energy = E(end);
    SUB_BEAM = defining_subbeam(BEAM_SORTED, par);
    div_x(end+1) = std(SUB_BEAM(:,end,4));
    div_y(end+1) = std(SUB_BEAM(:,end,5));
end


figure(11);
set(11, 'position', [50, 164, 600, 822]);
set(11, 'PaperPosition', [0., 0., 6, 8]);
set(11, 'color', 'w');
set(gca, 'fontsize', 16);
plot(E, div_x, E, div_y);
xlabel('E (GeV)  '), ylabel('Divergence (urad)  ');




