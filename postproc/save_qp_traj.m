function npart = save_qp_traj(istart, istop, n_3D_timestep_start, npt, datadir, fid)

my_SI_params;

myfile_rpinput = [datadir 'rpinput'];
n0 = my_get_quickpic_param(myfile_rpinput, 'Plasma_Density');
DT = my_get_quickpic_param(myfile_rpinput, 'DT');
DFPHA_BEAM = my_get_quickpic_param(myfile_rpinput, 'DFPHA_BEAM');
Box_X = my_get_quickpic_param(myfile_rpinput, 'Box_X');
Box_Y = my_get_quickpic_param(myfile_rpinput, 'Box_Y');

omega_p = 5.64e4*sqrt(n0);  % Plasma frequency in s-1

BEAMTAG_start = my_get_quickpic_beamtag(datadir, '01', sprintf('%.4d', n_3D_timestep_start));
alpha = max((BEAMTAG_start(:,2))); 
BEAMTAG_end = my_get_quickpic_beamtag(datadir, '01', sprintf('%.4d', n_3D_timestep_start+(npt-1)*DFPHA_BEAM));
REF_TAG_end = sort(alpha*BEAMTAG_end(:,1) + BEAMTAG_end(:,2));


BEAM_SORTED = zeros(istop-istart+1, npt, 6);

% npart_tmp = 0;
npart_tmp = size(REF_TAG_end, 1);
IND2 = 1:npart_tmp;
for i=1:npt
    n_3D_timestep_str = sprintf('%.4d', n_3D_timestep_start+(i-1)*DFPHA_BEAM);
    if mod(i,10)==0
        disp(n_3D_timestep_str);
    end
    BEAM = my_get_quickpic_phasespace(datadir, '01', n_3D_timestep_str); % auto-extracts beam part filenames exists
    
    BEAM(:,1) = BEAM(:,1) * SI_c / omega_p * 1e6; % to um
    BEAM(:,2) = BEAM(:,2) * SI_c / omega_p * 1e6; % to um
    BEAM(:,3) = BEAM(:,3) * SI_c / omega_p * 1e6; % to um
    BEAM(:,4) = BEAM(:,4) ./ BEAM(:,6) * 1e6; % to urad
    BEAM(:,5) = BEAM(:,5) ./ BEAM(:,6) * 1e6; % to urad
    BEAM(:,6) = BEAM(:,6) * SI_em * SI_c^2/SI_e / 1e9; % total energy [GeV]
    
    BEAMTAG = my_get_quickpic_beamtag(datadir, '01', n_3D_timestep_str); % auto-extracts beam part filenames exists
    NEWTAG = alpha*BEAMTAG(:,1) + BEAMTAG(:,2);  
    [SORTED_TAG, IND] = sort(NEWTAG);
    if size(NEWTAG,1) ~= npart_tmp
        [~,IND2,~] = intersect(SORTED_TAG, REF_TAG_end);
        npart_tmp = size(NEWTAG, 1);
    end
    BEAM = BEAM(IND(:),:);
    BEAM = BEAM(IND2(:),:);
%     [~,IND,~] = intersect(NEWTAG, REF_TAG_end);
%     BEAM = BEAM(IND(:),:);
    BEAM(:,3) = 1e6*SI_c*(i-1)*DFPHA_BEAM*DT/omega_p - BEAM(:,3);
    BEAM_SORTED(:,i,:) = BEAM(istart:istop, :);
end

cond2 = max(abs(BEAM_SORTED(:,:,1)),[],2) < 0.95*Box_X/2.;
cond3 = max(abs(BEAM_SORTED(:,:,2)),[],2) < 0.95*Box_Y/2.;
cond = cond2 & cond3;
BEAM_SORTED = BEAM_SORTED(cond,:,:);

disp('Number of macro-particles saved:');
disp(size(BEAM_SORTED,1));

fwrite(fid, permute(BEAM_SORTED, [3 2 1]), 'float');

npart = size(BEAM_SORTED,1);

end








