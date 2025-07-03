% Simulating B-mode Ultrasound Images Example
function PWNExMHz_ref_att0p1f2_BA6_400kPa_input(attRef,freqMHz,pressurekPa)
attRef = str2double(attRef);
freqMHz = str2double(freqMHz);
pressurekPa = str2double(pressurekPa);

addpath(genpath('/opt/MATLAB Add-Ons'))
parallel.gpu.enableCUDAForwardCompatibility(true)

% Execution
gpu          = true;
if gpu
    DATA_CAST       = 'gpuArray-single';
else
    DATA_CAST       = 'single';
end

% Properties of the propagation medium
c0          = 1540;                      % [m/s]
rho0        = 1000;                    % [kg/m^3]
alpha_power = 2;

% Properties of the transducer
element_pitch   = 0.3e-3;       % from verasonics specifications
element_length  = 4e-3;         % measured from L11-4v

% =========================================================================
%% DEFINE THE K-WAVE GRID
% =========================================================================
grid_size_x     = 5.5e-2;       % [m]
grid_size_y     = 4e-2;         % [m]
grid_size_z     = 4.5e-3;       % [m]
depth           = grid_size_x;  % imaging depth [m]
cfl             = 0.3;          % CFL number, could be 0.3 or 0.5
max_f0          = 18e6;

% calculate the grid spacing based on Nyquist and max_f0
dx = c0 / (2 * max_f0);   % [m]

% compute the size of the grid
Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dx);
Nz = roundEven(grid_size_z / dx);

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% create the time array
t_end           = depth*2/c0;     % [s];    % total compute time [s]
kgrid.makeTime(c0, cfl, t_end);

% =========================================================================
%% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
tone_burst_cycles = 10;

tone_burst_freq = freqMHz*1e6;
source_strength = pressurekPa*1e3;          % [Pa]
input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer (pitch = 0.3mm, length = 5 mm)
transducer.number_elements = 128;  	% total number of transducer elements
transducer.element_width = round(element_pitch/dx);             % width of each element [grid points]
transducer.element_length = round(element_length/dx);  	% length of each element [grid points]
transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = Inf;              % focus distance [m]
transducer.elevation_focus_distance = Inf;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = ones(transducer.number_elements, 1);
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);
transducer.properties;

% =========================================================================
%% DEFINE THE MEDIUM PROPERTIES
% =========================================================================
% Maps of UNIFORM B/A and AC
BonA_map = single(6*ones(Nx,Ny,Nz));
att_map = single(attRef/100*ones(Nx,Ny,Nz));

medium.sound_speed_ref = c0;
medium.sound_speed = c0;
medium.BonA = BonA_map;
medium.alpha_coeff = att_map;
medium.alpha_power = alpha_power;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================
% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', 'auto', ...
    'DataCast', DATA_CAST, 'PlotSim', false, 'SaveToDisk',false};

% loop through the scan lines
for aa = 1:2
    filename = fullfile(pwd,"random_media","BAPW_STD2_REF2025_"+aa+".mat");
    load(filename,'density');
    % update the command line status
    disp('');
    disp('Computing scan line PLANE WAVE...');

    % load the current section of the medium
    medium.density = density;
    if gpu
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});
    else
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, transducer, transducer, input_args{:});
    end
    fs = 1/kgrid.dt;
    rf_prebf(:,:,aa) = sensor_data';
end
file_out = ['PWNE',num2str(round(tone_burst_freq/1e6)), ...
    'MHz_ref_att0p',num2str(round(attRef),"%02d"),...
    'f',num2str(alpha_power*10),...
    '_BA6_nc', num2str(tone_burst_cycles),...
    '_',num2str(source_strength/1000),'kPa'];
disp(file_out)
pitch = transducer.element_width*dx;
save(file_out,'rf_prebf','fs','c0','pitch');

clear rf_prebf transducer

return