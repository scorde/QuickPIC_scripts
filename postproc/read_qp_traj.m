function BEAM_SORTED = read_qp_traj( datadir, istart, istop, npt )


fid = fopen([datadir 'qp_traj'], 'rb');
fseek(fid, (istart-1)*npt*6, 'bof');
npart = istop - istart + 1;
BEAM_SORTED = fread(fid, npart*npt*6, 'float');
fclose(fid);
BEAM_SORTED = reshape(BEAM_SORTED, 6, npt, npart);
BEAM_SORTED = permute(BEAM_SORTED, [3 2 1]);


end

