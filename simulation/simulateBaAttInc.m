function simulateBaAttInc(baseDir)
%% Simulating B-mode Ultrasound Images Example
addpath(genpath(pwd));
% set to 'single' or 'gpuArray-single' to speed up computations
DATA_CAST       = 'gpuArray-single';     
simuNames = {'BaAttInc1','BaAttInc2','BaAttInc3'};

for iSim = 1:length(simuNames)
% =========================================================================
%% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 40;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]

% set total number of grid points not including the PML
Nx = 1024 - 2 * pml_x_size;      % [grid points]
Ny = 675 - 2 * pml_y_size;      % [grid points]
Nz = 160 - 2 * pml_z_size;      % [grid points]

% calculate the spacing between the grid points
dx = 0.06e-3;
dy = dx;                        % [m]
dz = dx;                        % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
%% DEFINE THE MEDIUM
% =========================================================================
% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
stdDensity = 2/100;
medium.alpha_power = 2;

% Properties of the inclusion
radius_disk = (9)*1e-3;
center_depth = 22.5e-3;
switch iSim
    case 1
        BaBack = 6; BaInc = 9;
        alphaBack = 0.1; alphaInc = 0.1;

    case 2
        BaBack = 9; BaInc = 9;
        alphaBack = 0.1; alphaInc = 0.18;

    case 3
        BaBack = 6; BaInc = 9;
        alphaBack = 0.1; alphaInc = 0.18;
end

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% B over A map
medium.BonA = zeros(Nx,Ny,Nz);
for mm=1:Nz
    medium.BonA(:,:,mm) =(BaInc - BaBack) * makeDisc(Nx, Ny, ...
        round(center_depth/dx), Ny/2, round(radius_disk/dx));
end
medium.BonA = single(medium.BonA + BaBack);

% Attenuation map
medium.alpha_coeff = zeros(Nx,Ny,Nz);
for mm=1:Nz
    medium.alpha_coeff(:,:,mm) =(alphaInc - alphaBack) * makeDisc(Nx, Ny, ...
        round(center_depth/dx), Ny/2, round(radius_disk/dx));
end
medium.alpha_coeff = single(medium.alpha_coeff + alphaBack);

medium.sound_speed_ref = c0;
medium.sound_speed = c0;
medium.density = (ones(Nx,Ny,Nz) + stdDensity*randn(Nx,Ny,Nz))*rho0;

% =========================================================================
%% DEFINE THE INPUT SIGNAL
% =========================================================================
tone_burst_freq = 6e6;        % [Hz]
tone_burst_cycles = 10;
source_strength_vector = [80 400]*1e3;

% =========================================================================
%% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 128;  	% total number of transducer elements
transducer.element_width = 5;       % width of each element [grid points]
transducer.element_length = 100;  	% length of each element [grid points]
transducer.element_spacing = 0;  	% spacing between elements [grid points]
transducer.radius = inf;            % radius of curvature of transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, ...
    Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = Inf;              % [m]
transducer.elevation_focus_distance = Inf;    % [m]
transducer.steering_angle = 0;    % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = ones(transducer.number_elements, 1);

%% Looping two pressure levels
for ss = 1:length(source_strength_vector)
    %% RUN THE SIMULATION
    % append input signal used to drive the transducer
    source_strength =source_strength_vector(ss);          % [Pa]
    input_signal_norm = toneBurst(1/kgrid.dt, ...
        tone_burst_freq, tone_burst_cycles);
    input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

    % Transducer
    transducer.input_signal = input_signal;
    transducerSample = kWaveTransducer(kgrid, transducer);

    % RUN
    input_args = {'PMLInside', false,...
        'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
        'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, ...
        transducerSample, transducerSample, input_args{:});
    rf_prebf(:,:,ss) = sensor_data';
end

%% Saving
element_pitch = (transducer.element_width+transducer.element_spacing)*dy;
z = 0:size(rf_prebf,1)-1; z = z*kgrid.dt*c0/2;
x = 0:size(rf_prebf,2)-1; x = x-mean(x); x = x *element_pitch;
rf1 = rf_prebf(:,:,1);
rf2 = rf_prebf(:,:,2);
fs = 1/kgrid.dt;
save(fullfile(baseDir,['rf_prebf_',simuNames{iSim},'.mat']),...
    'rf1','rf2','x','z','fs','medium');
% focal_number = 3;
% rf1 = BFangle(rf1,transducer.number_elements,fs,c0, ...
%     element_pitch,'rect',focal_number,0);
% rf2 = BFangle(rf2,transducer.number_elements,fs,c0, ...
%     element_pitch,'rect',focal_number,0);
% save(fullfile(BaseDir,['rf_',simuNames{iSim},'.mat']),...
%     'rf1','rf2','x','z','fs','medium');

end

end
