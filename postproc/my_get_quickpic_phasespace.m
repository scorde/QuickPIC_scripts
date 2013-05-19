function qp_BEAM = my_get_quickpic_phasespace(mydir, n_beam_str, n_timestep_str)

qp_BEAM = [];

% establish quickpic output format
qp_version_suffix = my_get_quickpic_format(mydir);
if( strfind(qp_version_suffix, '.h5') )
  % newer QP versions: output in a single file per timestep
  myfile = [mydir 'RAW-BEAM/' num2str(n_beam_str, '%.2d') '/RAW-BEAM-' n_beam_str '_' n_timestep_str '.h5'];
  qp_BEAM = double(my_read_h5_beamdata(myfile));
else
% older QP versions: output split into multiple files per timestep

% writes the first line to get starting number (it seems to not always be 0 for beam 2)
temp_filename = [mydir 'PHA-BEAM/' n_beam_str '/' 'temp.txt'];
system(['ls -la ' mydir 'PHA-BEAM/' n_beam_str '/PHA-BEAM-' '*-' n_beam_str '_' n_timestep_str '* > ' temp_filename]);
% gets number of parts
fid = fopen(temp_filename, 'r');
ls_output = 0;
while( ls_output ~= -1 )
ls_output = fgets(fid);
  if( ls_output ~= -1)
    n_str = findstr('PHA-BEAM-', ls_output);
    if(n_str == [])
      warning(['EA: part file not found for ' 'PHA-BEAM-****' '-' n_beam_str '_' n_timestep_str '.hdf' '.  Stopping execution.']);
      stop;
    end% 
    part_filename = ls_output(n_str:n_str+24);
    myfile = [mydir 'PHA-BEAM/' n_beam_str '/'  part_filename];
    qp_BEAM_part = double(my_read_hdf(myfile));
    qp_BEAM = [qp_BEAM; qp_BEAM_part];
  end% if
end% while
fclose(fid);

end% of

return;
