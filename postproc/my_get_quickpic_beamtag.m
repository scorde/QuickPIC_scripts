function qp_tag = my_get_quickpic_beamtag(mydir, n_beam_str, n_timestep_str)

% establish quickpic output format
qp_version_suffix = my_get_quickpic_format(mydir);
if( strfind(qp_version_suffix, '.h5') )
  % newer QP versions: output in a single file per timestep
  myfile = [mydir 'RAW-BEAM/' num2str(n_beam_str, '%.2d') '/RAW-BEAM-' n_beam_str '_' n_timestep_str '.h5'];
  qp_tag = double(my_read_h5_beamtag(myfile));
end

return;