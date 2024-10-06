% Simulating B-mode Ultrasound Images Example


clearvars; close all; clc;
% simulation settings
% addpath(genpath('/opt/MATLAB Add-Ons'));
addpath(genpath(pwd));

DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 40;                % [grid points]
%pml_y_size = 10;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]

% set total number of grid points not including the PML
Nx = 1024 - 2 * pml_x_size;      % [grid points]
Ny = 675 - 2 * pml_y_size;      % [grid points]
Nz = 160 - 2 * pml_z_size;      % [grid points]

% set desired grid size in the x-direction not including the PML
%x = 60e-3;                      % [m]
%disp(['axial (mm): ',num2str(x*1e3),', lateral (mm): ',num2str(1e3*x*Ny/Nx),', elevation (mm): ',num2str(1e3*x*Nz/Nx)])

% calculate the spacing between the grid points
dx = 0.06e-3;
%dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
densityStd = 2/100;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);
medium.sound_speed_ref = c0;
medium.sound_speed = c0;

%%

% define properties of the input signal
for ii = 1:6 % number of sims
    %%
    % =========================================================================
    % DEFINE THE MEDIUM PROPERTIES
    % =========================================================================

    % define a large image size to move across
    %number_scan_lines = 500;
    Nx_tot = Nx;
    Ny_tot = Ny;
    Nz_tot = Nz;

    % Maps of B/A and AC
    switch ii
        case 1
            att = 0.4; ba = 6;
        case 2
            att = 0.4; ba = 9;
        case 3
            att = 0.4; ba = 12;
        case 4
            att = 0.6; ba = 6;
        case 5
            att = 0.6; ba = 9;
        case 6
            att = 0.6; ba = 12;
    end
    medium.alpha_coeff = single(att*ones(Nx,Ny_tot,Nz));
    medium.BonA = single(ba*ones(Nx,Ny_tot,Nz));
    medium.alpha_power = 1.2;
    medium.density = rho0*(1 + densityStd*randn(Nx,Ny_tot,Nz));

    % update the command line status
    disp('');
    disp('Computing scan line PLANE WAVE...');

    % load the current section of the medium

    param.c0 = c0;
    param.alpha_coeff = att;
    param.BonA = ba;
    param.alpha_power = medium.alpha_power;
    param.densityStd = densityStd;

    %%
    % =========================================================================
    % DEFINE THE INPUT SIGNAL
    % =========================================================================

    tone_burst_cycles = 3.5;
    tone_burst_freq = 5e6;        % [Hz]

    source_strength_vector = [80,240,400]*1e3;
    for ss = 1:length(source_strength_vector)
        % append input signal used to drive the transducer
        source_strength =source_strength_vector(ss);          % [Pa]
        input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq*1e6, tone_burst_cycles);
        input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

        % =========================================================================
        % DEFINE THE ULTRASOUND TRANSDUCER
        % =========================================================================

        % physical properties of the transducer
        transducer.number_elements = 128;  	% total number of transducer elements
        transducer.element_width = 5;       % width of each element [grid points]
        transducer.element_length = 100;  	% length of each element [grid points]
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
        %transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
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
        % RUN THE SIMULATION
        % =========================================================================
        % set the input settings
        [~,~,~] = mkdir(fullfile(pwd,'temp'));
        input_args = {...
            'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
            'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false,...
            'DataPath',fullfile(pwd,'temp')};
        
        for rr=1:4
            medium.density = rho0*(1 + densityStd*randn(Nx,Ny_tot,Nz));
            sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});
            rf_prebf(:,:,rr) = sensor_data';
        end
        fs = 1/kgrid.dt;

        file_out = ['PWNE_5MHz_homo_BA',...
            num2str(ba),'_att0p',num2str(att*10,'%.1i'),...
            'f1p2_nc3p5_',...
            num2str(source_strength/1000),'kPa'];
        save(file_out,'rf_prebf','fs','c0','param'); clear rf; clear transducer
    end
end

return

