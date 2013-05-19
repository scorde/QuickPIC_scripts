function WRITE_RP(myfilein, myfileout, input_struct)
% Generate rpinput file by replacing modified lines in a standard rpinput
% A bunch of nastiness, Don't touch!

fid = fopen(myfilein, 'r');
fidout = fopen(myfileout, 'w');

section_beam = 0;
section_plasma = 0;
section_neutral = 0;
section_species = 0;
section_simsys = 0;
section_nbeams = 0;
section_simtime = 0;
section_fielddiag = 0;
section_beamdiag = 0;
section_plasmadiag = 0;
section_beamphasediag = 0;
section_restart = 0;

row = 0;
n_row = 0;

while( row ~= -1 )

  n_row = n_row + 1;
  row = fgets(fid);
  
  if( row ~= -1 )
    % TAG SECTIONS
    if( ~isempty(strfind(row, '&Pipeline')) )
      section_pipeline = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_pipeline = 0;
    end% if
    if( ~isempty(strfind(row, '&Simulation_Sys')) )
      section_simsys = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_simsys = 0;
    end% if
    if( ~isempty(strfind(row, '&Num_Beams')) )
      section_nbeams = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_nbeams = 0;
    end% if
    if( ~isempty(strfind(row, '&Beam')) && isempty(strfind(row, '&Beam_')) )
      section_beam = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_beam = 0;
    end% if
    if( ~isempty(strfind(row, '&Plasma')) && isempty(strfind(row, '&Plasma_')) )
      section_plasma = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_plasma = 0;
    end% if
    if( ~isempty(strfind(row, '&Neutral')) && isempty(strfind(row, '&Neutral_')) )
      section_neutral = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_neutral = 0;
    end% if
    if( ~isempty(strfind(row, '&Species')) )
      section_species = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_species = 0;
    end% if
    if( ~isempty(strfind(row, '&Simulation_time')) )
      section_simtime = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_simtime = 0;
    end% if
    if( ~isempty(strfind(row, '&Field_Diag')) )
      section_fielddiag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_fielddiag = 0;
    end% if
    if( ~isempty(strfind(row, '&Beam_Diag')) )
      section_beamdiag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_beamdiag = 0;
    end% if
    if( ~isempty(strfind(row, '&Plasma_Diag')) )
      section_plasmadiag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_plasmadiag = 0;
    end% if
    if( ~isempty(strfind(row, '&Beam_Phase_Space_Diag')) )
      section_beamphasediag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_beamphasediag = 0;
    end% if
    if( ~isempty(strfind(row, '&Restart_File')) )
      section_restart = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_restart = 0;
    end% if
    

    % STAGES SECTION
    if( ~isempty(strfind(row, 'Num_Stages') ) && section_pipeline == 1)
       fprintf(fidout, [' Num_Stages = ' num2str(input_struct.comp.num_stages) '\n']  );

    % SIMSYS SECTION
    elseif( ~isempty(strfind(row, 'Box_X') ) && section_simsys == 1)
       fprintf(fidout, [' Box_X=' num2str(input_struct.size.Box_X) ', Box_Y=' num2str(input_struct.size.Box_Y) ', Box_Z=' num2str(input_struct.size.Box_Z) ',' '\n']  );
    elseif( ~isempty(strfind(row, 'INDX') ) && section_simsys == 1)
       fprintf(fidout, [' INDX =  ' num2str(input_struct.size.INDX) ' , INDY = ' num2str(input_struct.size.INDY) ', INDZ = ' num2str(input_struct.size.INDZ) '\n']  );

    % NUMBEAMS SECTION
    elseif( ~isempty(strfind(row, 'NBeams') ) && section_nbeams == 1)
       fprintf(fidout, [' NBeams = ' num2str(1) '\n']  );

    % BEAM SECTION(S)
    elseif( section_beam == 1)
       %for(n_beam=1:N_beams),
          fprintf(fidout, '\n'  );
          fprintf(fidout, ['&Beam' '\n']  );
          fprintf(fidout, [' BEAM_EVOLUTION = .true.' '\n']  );
          fprintf(fidout, [' MIN_BEAM_PARTICLE = 8' '\n']  );
          fprintf(fidout, [' NPX =  ' num2str(2^7) ', NPY = ' num2str(2^7) ', NPZ = ' num2str(2^8) '\n']);
          fprintf(fidout, [' Charge = ' num2str(input_struct.beam.charge, '%.1f') '\n']  );
          qp_par_mass_str = num2str(input_struct.beam.mass, '%.3G');
          qp_par_mass_str = strrep(qp_par_mass_str,'E+0','E');
          qp_par_mass_str = strrep(qp_par_mass_str,'E+','E');
          fprintf(fidout, [' Mass = ' qp_par_mass_str '\n']  );
          qp_gam_part_str = num2str(input_struct.beam.gamma, '%.3G');
          qp_gam_part_str = strrep(qp_gam_part_str,'E+0','E');
          qp_gam_part_str = strrep(qp_gam_part_str,'E+','E');
          fprintf(fidout, [' Gamma = ' qp_gam_part_str ',' '\n']  );
          qp_num_part_str = num2str(input_struct.beam.N_particles, '%.3G');
          qp_num_part_str = strrep(qp_num_part_str,'E+0','E');
          qp_num_part_str = strrep(qp_num_part_str,'E+','E');
          fprintf(fidout, [' Num_Particle = ' qp_num_part_str ',' '\n']);
          fprintf(fidout, [' VDX =   0.0, VDY =   0.0, VDZ =  0.0' '\n']  );
          fprintf(fidout, [' Init_Routine = 1' '\n']  );
          fprintf(fidout, [' BEAM_PROFILE = ''test.hdf''' '\n']  );
          fprintf(fidout, [' QUIET_START = .true.' '\n']  );
          fprintf(fidout, [' Parameter_Array(1:1,1:3) = ' num2str(input_struct.pos.X_center, '%.2f') ',' num2str(input_struct.pos.Y_center, '%.2f')  ',' num2str(input_struct.pos.Z_center, '%.2f') '\n']);
          fprintf(fidout, [' Parameter_Array(2:2,1:3) = ' num2str(input_struct.beam.sigma_x, '%.2f') ',' num2str(input_struct.beam.sigma_y, '%.2f')  ',' num2str(input_struct.beam.sigma_z, '%.2f') '\n']);
          fprintf(fidout, [' Parameter_Array(3:3,1:3) = ' num2str(input_struct.beam.emit_x, '%.3f') ',' num2str(input_struct.beam.emit_y, '%.3f') ',' num2str(0., '%.3f')  '\n'] );
          fprintf(fidout, [' Parameter_Array(4:4,1:3) = ' num2str(0., '%.3f') ',' num2str(-0., '%.3f') ',' num2str(0., '%.3f')  '\n'] );
          fprintf(fidout, [' Parameter_Array(5:5,1:3) = ' num2str(0., '%.3f') ',' num2str(-0., '%.3f') ',' num2str(0., '%.3f')  '\n'] );
          fprintf(fidout, [' Use_Shifter = .false.' '\n']  );
          fprintf(fidout, [' Shifter_Nsec = 4' '\n']  );
          fprintf(fidout, [' Shifter_Parameter(1:1,1:4) = 0.,0.,1.5,0.' '\n']  );
          fprintf(fidout, [' Shifter_Parameter(2:2,1:4) = 0.,0.,0.,0.' '\n']  );
          fprintf(fidout, [' Shifter_Parameter(3:3,1:4) = 0.,78.,155.1,155.2' '\n']  );
          fprintf(fidout, [' Use_Destroyer = .true.' '\n']  );
          fprintf(fidout, [' Destroyer_NCriteria = 5' '\n']  );
          fprintf(fidout, [' Destroyer_Criteria(1:1,1:5)=1,1,2,2,6' '\n']  );
          fprintf(fidout, [' Destroyer_Criteria(2:2,1:5)=0,' num2str(input_struct.size.Box_X-2, '%.0f') ',0,' num2str(input_struct.size.Box_Y-2, '%.0f')  ',0' '\n']);
          fprintf(fidout, [' Destroyer_Criteria(3:3,1:5)=2,' num2str(input_struct.size.Box_X, '%.0f') ',2,' num2str(input_struct.size.Box_Y, '%.0f')  ',100' '\n']);
          fprintf(fidout, [' Use_Radiation_Damping = .false.' '\n']  );
          fprintf(fidout, ['/' '\n']  );
      %end% for
      section_beam = 0;


      % PLASMA SECTION
      elseif( ~isempty(strfind(row, 'Plasma_Density') ) && section_plasma == 1)
         Plasma_Density_str = num2str(input_struct.plasma.density, '%.3G');
         Plasma_Density_str = strrep(Plasma_Density_str,'E+0','E');
         Plasma_Density_str = strrep(Plasma_Density_str,'E+','E');
         fprintf(fidout, [' Plasma_Density=' Plasma_Density_str '\n']  );
      elseif( ~isempty(strfind(row, 'Nspecies') ) && section_plasma == 1)
         if( input_struct.plasma.PREION == 1)
            fprintf(fidout, [' Nspecies=1' '\n']  );
         else
            fprintf(fidout, [' Nspecies=0' '\n']  );
         end% if
      elseif( ~isempty(strfind(row, 'Nneutral') ) && section_plasma == 1)
         if( input_struct.plasma.PREION == 1)
            fprintf(fidout, [' Nneutrals=0' '\n']  );
         else
            fprintf(fidout, [' Nneutrals=1' '\n']  );
         end% if

      % SPECIES SECTION 
      elseif( ~isempty(strfind(row, 'NP2') ) && section_species == 1)
         n_plasma_particles_per_cell = 4; % form UCLA
         fprintf(fidout, [' NP2 = ' num2str(round(sqrt((2^input_struct.size.INDX)^2 * n_plasma_particles_per_cell))) '\n']  );
      elseif( ~isempty(strfind(row, 'Charge') ) && section_species == 1)
         fprintf(fidout, [' Charge = ' num2str( input_struct.plasma.charge, '%.1f') '\n'] );
      elseif( ~isempty(strfind(row, 'Mass') ) && section_species == 1)
         fprintf(fidout, [' Mass = ' num2str( input_struct.plasma.mass, '%.1f') '\n'] );
      elseif( ~isempty(strfind(row, 'Profile_type') ) && section_species == 1)
          if input_struct.plasma.profile == 0
              fprintf(fidout, ' Profile_type=0 \n' );
          elseif input_struct.plasma.profile == 1
              fprintf(fidout, ' Profile_type=21 \n' );
          end
      elseif( ~isempty(strfind(row, 'Prof_Nsec') ) && section_species == 1)
          fprintf(fidout, [' Prof_Nsec = ' num2str( input_struct.plasma.n_point) '\n'] );
      elseif( ~isempty(strfind(row, 'Prof_Parameter(1,') ) && section_species == 1)
          fprintf(fidout, [' Prof_Parameter(1,1:' num2str(input_struct.plasma.n_point) ') = '] );
          for i=1:(input_struct.plasma.n_point-1)
              fprintf(fidout, [num2str(input_struct.plasma.n(i),'%.2f') ',']);
          end
          fprintf(fidout,[num2str(input_struct.plasma.n(input_struct.plasma.n_point),'%.2f') '\n']);
      elseif( ~isempty(strfind(row, 'Prof_Parameter(2,') ) && section_species == 1)
          fprintf(fidout, [' Prof_Parameter(2,1:' num2str(input_struct.plasma.n_point) ') = '] );
          for i=1:(input_struct.plasma.n_point-1)
              fprintf(fidout, [num2str(input_struct.plasma.r(i),'%.2f') ',']);
          end
          fprintf(fidout,[num2str(input_struct.plasma.r(input_struct.plasma.n_point),'%.2f') '\n']);

      % NEUTRAL SECTION
      elseif( ~isempty(strfind(row, 'Neutral_gas') ) && section_neutral == 1)
         fprintf(fidout, [' Neutral_gas = ' num2str(input_struct.plasma.Z) '\n']  );
         
      % SIMTIME SECTION
      elseif( ~isempty(strfind(row, 'TEND') ) && section_simtime == 1) 
         TEND_str = num2str(input_struct.time.TEND, '%10G');
         TEND_str = strrep(TEND_str,'E+0','E');
         TEND_str = strrep(TEND_str,'E+','E');
         fprintf(fidout, [' TEND =' TEND_str ', DT = ' num2str(input_struct.time.DT, '%.1f') '  ,' '\n']  );

      % DIAG SECTION
      elseif( ~isempty(strfind(row, 'DFESLICE') ) && section_fielddiag == 1) 
         fprintf(fidout, [' DFESLICE=' num2str(input_struct.time.DT_OUTPUT) ', EX0=0, EY0=' num2str(input_struct.pos.Y_center, '%.0f') ', EZ0=0' '\n']  );
      elseif( ~isempty(strfind(row, 'DFBSLICE') ) && section_fielddiag == 1) 
         fprintf(fidout, [' DFBSLICE=' num2str(input_struct.time.DT_OUTPUT) ', BX0=0, BY0=' num2str(input_struct.pos.Y_center, '%.0f') ', BZ0=0' '\n']  );
      elseif( ~isempty(strfind(row, 'DFQEB')) && isempty(strfind(row, 'DFQEBSLICE'))  && section_beamdiag == 1)
         if(input_struct.diag.store_QEB_3D)
            fprintf(fidout, [' DFQEB=' num2str(input_struct.time.DT_OUTPUT) ',' '\n']  );
         else
            fprintf(fidout, [' DFQEB=0,' '\n']  );
         end% if
      elseif( ~isempty(strfind(row, 'DFQEBSLICE'))  && section_beamdiag == 1)
         fprintf(fidout, [' DFQEBSLICE=' num2str(input_struct.time.DT_OUTPUT) ' , QEBX0=0., QEBY0=' num2str(input_struct.pos.Y_center, '%.0f') ',  QEBZ0=0' '\n']  );
      elseif( ~isempty(strfind(row, 'DFQEPSLICE'))  && section_plasmadiag == 1)
         fprintf(fidout, [' DFQEPSLICE=' num2str(input_struct.time.DT_OUTPUT) ' , QEPX0= 0, QEPY0=' num2str(input_struct.pos.Y_center, '%.0f') ',  QEPZ0=0' '\n']  );
      elseif( ~isempty(strfind(row, 'DUMP_PHA_BEAM'))  && section_beamphasediag == 1)
         fprintf(fidout, [' DUMP_PHA_BEAM=.true., DFPHA_BEAM=' num2str(input_struct.time.DT_OUTPUT) ',' '\n']  );
      elseif( ~isempty(strfind(row, 'DSAMPLE_BEAM'))  && section_beamphasediag == 1)
         fprintf(fidout, [' DSAMPLE_BEAM = ' num2str(2^5) '\n']  );

      % RESTART SECTION
      elseif( ~isempty(strfind(row, 'DUMP_RST_FILE'))  && section_restart == 1)
         if( input_struct.time.DT_RESTART)
            fprintf(fidout, [' DUMP_RST_FILE = .true.,  DFRST=' num2str(input_struct.time.DT_RESTART) '\n']  );
         else
            fprintf(fidout, [' DUMP_RST_FILE = .false.,  DFRST=1' '\n']  );
         end% if
      else
         fprintf(fidout, row);
    end%if
  end%if
end% while
fclose(fid);
fclose(fidout);
