function param_struct = CALC_RP(input_struct)

% import standard SI constants
eval(['run ' pwd '/SI_consts.m']);



%%%%%%%%%%%%%%%%%%%%%
% PLASMA ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%

param_struct.plasma.density  = input_struct.plasma.density; % plasma density [cm^-3]
param_struct.plasma.charge   = input_struct.plasma.charge;  % plasma particle charge [e]
param_struct.plasma.mass     = input_struct.plasma.mass;    % plasma particle mass [e_m]
param_struct.plasma.PREION   = input_struct.plasma.PREION;  % 0: non-ionized plasma, 1: pre-ionized plasma
param_struct.plasma.Z        = input_struct.plasma.Z;       % ion number
param_struct.plasma.profile  = input_struct.plasma.profile; % 0: uniform plasma, 1: hollow channel plasma

% Calc plasma parameters
param_struct.plasma.omega_p  = sqrt(param_struct.plasma.density*...
    1e6*SI_e^2/(SI_em*SI_eps0));                                                    % plasma frequency  [rad/s]
param_struct.plasma.lambda_p = 2*pi*SI_c*1e6/param_struct.plasma.omega_p;           % plasma wavelength [um]
param_struct.plasma.k_p      = param_struct.plasma.omega_p/SI_c;                    % plasma wavenumber [1/m]
param_struct.plasma.SD       = 1e6*SI_c/param_struct.plasma.omega_p;                % plasma skin depth [um]
param_struct.plasma.time     = 1/param_struct.plasma.omega_p;                       % characteristic time scale [s]
param_struct.plasma.rqm      = param_struct.plasma.mass/param_struct.plasma.charge; % mass to charge ratio
param_struct.plasma.field    = SI_em*SI_c*param_struct.plasma.omega_p/(1e9*SI_e);   % GV/m



%%%%%%%%%%%%%%%%%%%
% BEAM ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%

param_struct.beam.charge      = input_struct.beam.charge;                        % beam particle charge [e]
param_struct.beam.mass        = input_struct.beam.mass;                          % beam particle mass [e_m]
param_struct.beam.N_particles = input_struct.beam.N_particles;                   % number of beam particles
param_struct.beam.rqm         = param_struct.beam.mass/param_struct.beam.charge; % mass to charge ratio

% Calc gamma if not specified
if input_struct.beam.gamma == 0
   if input_struct.beam.energy == 0
      error('Must specify either gamma or energy');
   end
   param_struct.beam.gamma  = input_struct.beam.energy/...
       (input_struct.beam.mass*SI_eM/1e3);                  % beam gamma
   param_struct.beam.energy = input_struct.beam.energy;     % beam energy [GeV]
end

% Calc energy if not specified
if input_struct.beam.energy == 0
   if input_struct.beam.gamma == 0
      error('Must specify either gamma or energy');
   end
   param_struct.beam.energy = input_struct.beam.gamma*...
       input_struct.beam.mass*SI_eM/1e3;                    % beam energy [GeV]
   param_struct.beam.gamma  = input_struct.beam.gamma;      % beam gamma
end

% Betatron frequecy, wavelength
param_struct.beam.omega_b  = sqrt(param_struct.plasma.density*1e6*SI_e^2/... 
    (2*param_struct.beam.gamma*param_struct.beam.mass*SI_em*SI_eps0));     % betatron frequency [rad/s]
param_struct.beam.lambda_b = 2*pi*SI_c*1e6/param_struct.beam.omega_b;        % betatron wavelength [um]

% Calc beam size if beam_match == 1
if( input_struct.beam.beam_match )
  param_struct.beam.sigma_x = 1e3*sqrt(input_struct.beam.emit_x/...
      param_struct.plamsa.k_p*sqrt(2/param_struct.beam.gamma));      % beam size X [um]
  
  param_struct.beam.sigma_y = 1e3*sqrt(input_struct.beam.emit_y/...
      param_struct.plasma.k_p*sqrt(2/param_struct.beam.gamma));      % beam size Y [um]
  
  param_struct.beam.emit_x  = input_struct.beam.emit_x;              % beam norm X emmitance [mm*mrad]
  param_struct.beam.emit_y  = input_struct.beam.emit_x;              % beam norm Y emmitance [mm*mrad]

% Calc emittance if emit_match == 1
elseif( input_struct.beam.emit_match )
  param_struct.beam.emit_x  = 1e-6*param_struct.plasma.k_p*...
      (input_struct.beam.sigma_x)^2*sqrt(param_struct.beam.gamma/2); % beam norm X emmitance [mm*mrad]
  param_struct.beam.emit_y  = 1e-6*param_struct.plasma.k_p*...
      (input_struct.beam.sigma_y)^2*sqrt(param_struct.beam.gamma/2); % beam norm Y emmitance [mm*mrad]
  
  param_struct.beam.sigma_x = input_struct.beam.sigma_x;             % beam size X [um]
  param_struct.beam.sigma_y = input_struct.beam.sigma_y;             % beam size Y [um]

else
  param_struct.beam.sigma_x = input_struct.beam.sigma_x;             % beam size X [um]
  param_struct.beam.sigma_y = input_struct.beam.sigma_y;             % beam size Y [um]
  
  param_struct.beam.emit_x  = input_struct.beam.emit_x;              % beam norm X emmitance [mm*mrad]
  param_struct.beam.emit_y  = input_struct.beam.emit_x;              % beam norm Y emmitance [mm*mrad]
end

% Normalized beam divergence
param_struct.beam.angle_x = param_struct.beam.emit_x/param_struct.beam.sigma_x; % emit_x/sigma_x
param_struct.beam.angle_y = param_struct.beam.emit_y/param_struct.beam.sigma_y; % emit_x/sigma_x

% Calc sigma_z if z_match == 1
if( input_struct.beam.z_match )
   param_struct.beam.sigma_z = 1e6*sqrt(2)*param_struct.plasma.SD;   % beam size Z [um]
else
   param_struct.beam.sigma_z = input_struct.beam.sigma_z;            % beam size Z [um]
end

% Calc beam density
param_struct.beam.density = (10000^3)*param_struct.beam.N_particles/...
    ((2*pi)^(3/2)*param_struct.beam.sigma_x*param_struct.beam.sigma_y*...
    param_struct.beam.sigma_z);                                           % beam density [cm^-3]

% Calc N_b/N_p
param_struct.beam.ratio = param_struct.beam.density/param_struct.plasma.density; % beam density [n0]

% Calc k_p*sigma_x,y,z
param_struct.beam.kpsx = param_struct.beam.sigma_x/param_struct.plasma.SD; % beam length X [skin depths]
param_struct.beam.kpsy = param_struct.beam.sigma_y/param_struct.plasma.SD; % beam length Y [skin depths]
param_struct.beam.kpsz = param_struct.beam.sigma_z/param_struct.plasma.SD; % beam length Z [skin depths]

% Calc peak gaussian current
param_struct.beam.I_peak = 1e3*param_struct.beam.N_particles*...
    SI_e*SI_c/(sqrt(2*pi)*param_struct.beam.sigma_z);            % beam current [kA]

% Calc max bubble radius
param_struct.plasma.R_bubble = (1/0.84)*2*param_struct.plasma.SD*...
    sqrt(param_struct.beam.I_peak/17);                               % plasma bubble radius[um]



%%%%%%%%%%%%%%%%%%%
% SIZE ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%

% Determine simulation box size
Box_X_exact = max(input_struct.size.X_bubbles*param_struct.plasma.R_bubble,...
    input_struct.size.X_bunches*param_struct.beam.sigma_x);                     % [um]
Box_Y_exact = max(input_struct.size.X_bubbles*param_struct.plasma.R_bubble,...
    input_struct.size.X_bunches*param_struct.beam.sigma_x);                     % [um]
Box_Z_exact = max(input_struct.size.Z_waves*param_struct.plasma.lambda_p,...
    input_struct.size.Z_bunches*param_struct.beam.sigma_z);                    % [um]

% Round to nearet micron
param_struct.size.Box_X = ceil(Box_X_exact/10)*10; % [um]
param_struct.size.Box_Y = ceil(Box_Y_exact/10)*10; % [um]
param_struct.size.Box_Z = ceil(Box_Z_exact/10)*10; % [um]

% Determine grid spacing in terms of skin depth
skin_frac = 1/20;

% Determine number of cells, must be at least 2^6 per dim
param_struct.size.INDX = max(6,floor(log(param_struct.size.Box_X*1e-6*...
    param_struct.plasma.k_p/skin_frac)/log(2)));
param_struct.size.INDY = max(6,floor(log(param_struct.size.Box_Y*1e-6*...
    param_struct.plasma.k_p/skin_frac)/log(2)));
param_struct.size.INDZ = max(6,floor(log(param_struct.size.Box_Z*1e-6*...
    param_struct.plasma.k_p/skin_frac)/log(2)));

% Increase granularity if desired
param_struct.size.INDX = param_struct.size.INDX+input_struct.size.x_grain;
param_struct.size.INDY = param_struct.size.INDY+input_struct.size.y_grain;
param_struct.size.INDZ = param_struct.size.INDZ+input_struct.size.z_grain;

% Determine cell size
param_struct.size.Cell_X = param_struct.size.Box_X/(2^(param_struct.size.INDX)); % [um]
param_struct.size.Cell_Y = param_struct.size.Box_Y/(2^(param_struct.size.INDY)); % [um]
param_struct.size.Cell_Z = param_struct.size.Box_Z/(2^(param_struct.size.INDZ)); % [um]

% Do not allow transverse cell size to be smaller than longitudinal cell size
while param_struct.size.Cell_Z > param_struct.size.Cell_X
   param_struct.size.INDZ = param_struct.size.INDZ + 1;
   param_struct.size.Cell_Z = param_struct.size.Box_Z/(2^(param_struct.size.INDZ));
end

% Determine number of grid points
param_struct.size.Grid_X = 2^param_struct.size.INDX;
param_struct.size.Grid_Y = 2^param_struct.size.INDY;
param_struct.size.Grid_Z = 2^param_struct.size.INDZ;

% Cell Fraction
param_struct.size.Frac_X = param_struct.size.Cell_X/param_struct.plasma.SD;
param_struct.size.Frac_Y = param_struct.size.Cell_Y/param_struct.plasma.SD;
param_struct.size.Frac_Z = param_struct.size.Cell_Z/param_struct.plasma.SD;



%%%%%%%%%%%%%%%%%%
% PLASMA PROFILE %
%%%%%%%%%%%%%%%%%%

if param_struct.plasma.profile == 1    
    param_struct.plasma.n_point = input_struct.plasma.n_point; % number of points used to define channel
    param_struct.plasma.radius  = input_struct.plasma.radius;  % channel radius
    param_struct.plasma.width   = input_struct.plasma.width;   % annulus width
    
    % Calc plasma profile
    param_struct.plasma.r = linspace(0,floor(param_struct.size.Box_X/2),param_struct.plasma.n_point); % N radial points
    %param_struct.plasma.n = param_struct.plasma.density * ...
        %exp(-(param_struct.plasma.r - param_struct.plasma.radius).^2/(2*param_struct.plasma.width^2)); % radial density
    param_struct.plasma.n = exp(-(param_struct.plasma.r - param_struct.plasma.radius).^2/(2*param_struct.plasma.width^2));
    %param_struct.plasma.n = param_struct.plasma.n .* (param_struct.plasma.n > 1);
end

%%%%%%%%%%%%%%%%%%%
% TIME ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%

% Determine slice time
if input_struct.sim.DT == 0
    param_struct.time.DT = round(sqrt(2*param_struct.beam.gamma)/10 / 1.5); % 1.5 is deceleration factor
else
    param_struct.time.DT = input_struct.sim.DT;
end

if( input_struct.sim.BEAM_EV )
  param_struct.time.TEND       = floor(input_struct.sim.prop / (SI_c / param_struct.plasma.omega_p))+0.1;
  param_struct.time.DT_OUTPUT  = input_struct.sim.dump_freq; % output data every n'th timestep
  param_struct.time.DT_RESTART = param_struct.time.DT_OUTPUT*4;
else
  param_struct.time.TEND       = param_struct.time.DT+0.1;
  param_struct.time.DT_OUTPUT  = 1; % output data every timestep
  param_struct.time.DT_RESTART = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%
% POSITION ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%%%

% Determine beam position
param_struct.pos.X_center = round(param_struct.size.Box_X*input_struct.size.X_center);
param_struct.pos.Y_center = round(param_struct.size.Box_Y*input_struct.size.Y_center);
param_struct.pos.Z_center = round(param_struct.size.Box_Z*input_struct.size.Z_center);



%%%%%%%%%%%%%%%%%%%%%%%%%
% DIAGNOSTIC ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%%%%%

param_struct.diag.store_QEB_3D = input_struct.diag.store_QEB_3D;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONAL ATTRIBUTES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if input_struct.sim.BEAM_EV == 0
    param_struct.comp.num_stages = 1;
    param_struct.comp.run_time   = 1;
    param_struct.comp.mem        = 1024;
    param_struct.comp.tasks      = 8;
elseif input_struct.sim.BEAM_EV == 1
    param_struct.comp.num_stages = 2;
    param_struct.comp.run_time   = input_struct.sim.run_time;
    param_struct.comp.mem        = 4096;
    param_struct.comp.tasks      = 128;
end