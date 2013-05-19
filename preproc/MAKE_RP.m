% QuickPIC Matlab rpinput generation example script
% S. Gessner Sep 07, 2012

%clear all;

% import standard SI constants
SI_consts;

% specify template
rpinput_template_file = [pwd '/rpinputs/rpinput_template'];

% specify output
date_dir = GET_DATE_DIR;

rpinput_dir = [pwd '/rpinputs/' date_dir];
if ~exist(rpinput_dir,'dir')
    mkdir(rpinput_dir);
end

param_dir = [pwd '/params/' date_dir];
if ~exist(param_dir,'dir')
    mkdir(param_dir);
end

command_dir = [pwd '/commands/' date_dir];
if ~exist(command_dir,'dir')
    mkdir(command_dir);
end

rpinput_output_name = 'profile';
rpinput_output_file = [rpinput_dir 'rpinput_' rpinput_output_name];

write = 1;
% check to see if you want to overwrite file
if exist(rpinput_output_file,'file')
   reply = input(['File ' rpinput_output_file ' exists. \n Do you want to overwrite? y/n '], 's');
   if (strcmp(reply,'n'))
      disp('Ok. That''s cool.');
      write = 0;
   end
end


% INPUT TO RPINPUT

% simulation parameters
input_struct.sim.BEAM_EV       = 1;           % 0 : calc wake only, 1 : propagate and evolve beam
input_struct.sim.prop          = 0.1;         % propagation length of the beam [m]
input_struct.sim.DT            = 0;           % Delta T between beam pushes [1/omega_p]. If 0: use calc from formula
input_struct.sim.dump_freq     = 1;           % Dump frequency
input_struct.sim.run_time      = 1;           % Amount of computer time to run sim for, 1 if BEAM_EV = 0

% plasma parameters
input_struct.plasma.density    = 6e16;        % /cm^3
input_struct.plasma.charge     = -1.0;        % -1 for electron, +1 for positron
input_struct.plasma.mass       = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.plasma.PREION     = 1;           % 0 : non-ionized plasma 1: pre-ionized plasma
input_struct.plasma.Z          = 3;           % atomic number of plasma gas
input_struct.plasma.profile    = 1;           % 0: uniform plasma, 1: hollow channel plasma
input_struct.plasma.n_point    = 100;         % number of points used to create plasma profile
input_struct.plasma.radius     = 25;          % channel radius [um]
input_struct.plasma.width      = 3;           % annulus width [um]

% beam parameters
input_struct.beam.charge       = +1.0;        % -1 for electron, +1 for positron
input_struct.beam.mass         = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.beam.N_particles  = 1.0e10;      % Number of beam particles
input_struct.beam.gamma        = 40000;       % relativistic factor gamma, if 0 energy specified below
input_struct.beam.energy       = 0;           % beam mean energy [GeV], if 0 use gamma to calculate energy
input_struct.beam.sigma_x      = 10;          % Gaussian sigma_x [um]
input_struct.beam.sigma_y      = 10;          % Gaussian sigma_y [um]
input_struct.beam.sigma_z      = 14;          % Gaussian sigma_z [um]
input_struct.beam.emit_x       = 5.00;        % normalized X emittance [mm*mrad i.e. 1e-6]
input_struct.beam.emit_y       = 5.00;        % normalized Y emittance [mm*mrad i.e. 1e-6]
input_struct.beam.beam_match   = 0;           % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
input_struct.beam.emit_match   = 0;           % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
input_struct.beam.z_match      = 0;           % 1: override sigma_z with sqrt(2)/k_p, 0: do nothing

% size parameters
input_struct.size.Z_waves      = 2.;          % set box length by number of plasma wavelengths
input_struct.size.Z_bunches    = 0;           % set box length by bunch lengths
input_struct.size.X_bubbles    = 5;           % set box width by number of bubble radii
input_struct.size.X_bunches    = 0.0;         % set box width by number of bunch radii
input_struct.size.X_center     = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.size.Y_center     = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.size.Z_center     = 0.20;        % place bunch as fraction of box length, 0 at start of box, 1 at end
input_struct.size.x_grain      = 0;           % increase granularity in the x dimension by 2^x_grain
input_struct.size.y_grain      = 0;           % increase granularity in the y dimension by 2^y_grain
input_struct.size.z_grain      = 0;           % increase granularity in the z dimension by 2^z_grain

% diagnostic parameters
input_struct.diag.store_QEB_3D = 0;           % store full 3D beam phase space?

% RPINPUT CALCULATOR
param_struct = CALC_RP(input_struct);

% RPINPUT WRITER
if write
    WRITE_RP(rpinput_template_file, rpinput_output_file, param_struct);
    save([param_dir 'param_' rpinput_output_name '.mat'], 'param_struct');
    run_dir = WRITE_CMD(command_dir, rpinput_output_name, param_struct.comp.mem,...
        param_struct.comp.tasks, param_struct.comp.run_time);
end

%exit;
