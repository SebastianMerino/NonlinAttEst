% Simulating B-mode Ultrasound Images Example


clearvars; close all; clc;
addpath(genpath('/opt/MATLAB Add-Ons'))
% simulation settings
DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
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

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
tone_burst_freq_vector = [5e6];        % [Hz]
source_strength_vector = [400]*1e3;
alpha_power = 1.2;

for attInc = [10,14]
    for iFreq = 1:length(tone_burst_freq_vector)
        tone_burst_freq = tone_burst_freq_vector(iFreq);
        tone_burst_cycles = 10;
        for ss = 1:length(source_strength_vector)
            % append input signal used to drive the transducer
            source_strength =source_strength_vector(ss);          % [Pa]
            input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
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
            % DEFINE THE MEDIUM PROPERTIES
            % =========================================================================

            % define a large image size to move across
            Nx_tot = Nx;
            Ny_tot = Ny;
            Nz_tot = Nz;

            %% Maps of NONUNIFORM B/A and AC
            BonA_map = single(6*ones(Nx,Ny_tot,Nz));
            att_map = single(0.1*ones(Nx,Ny_tot,Nz));
            radius_disk = (9)*1e-3;
            center_depth = 22.5e-3;
            BonAdiff = 6;
            attDiff = (attInc/100 - 0.1);
            for mm=1:Nz
                BonA_map(:,:,mm) =single(BonAdiff * makeDisc(Nx_tot, Ny_tot, round(center_depth/dx), Ny_tot/2, round(radius_disk/dx)));
                att_map(:,:,mm) =single(attDiff * makeDisc(Nx_tot, Ny_tot, round(center_depth/dx), Ny_tot/2, round(radius_disk/dx)));
            end
            BonA_map = BonA_map + 6;
            att_map = att_map + 0.1;


            medium.sound_speed_ref = 1540;
            medium.sound_speed = 1540;

            % =========================================================================
            % RUN THE SIMULATION
            % =========================================================================
            % set the input settings
            input_args = {...
                'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
                'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

            % loop through the scan lines
            for aa = 1:4
                % for scan_line_index = 1% 1:number_scan_lines
                filename = [pwd,'/random_media/BAPW_STD2_SAM2023_',num2str(aa)];
                load(filename);
                % update the command line status
                disp('');
                disp(['Computing scan line PLANE WAVE...']);

                % load the current section of the medium
                medium.BonA = BonA_map;
                medium.alpha_coeff = att_map;
                medium.alpha_power = alpha_power;
                medium.density = density;

                param.c0 = c0;
                param.att_map_setting = att_map(:,:,1);
                param.BonA_map_setting = BonA_map(:,:,1);
                param.att_power_law = medium.alpha_power;

                sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});
                fs = 1/kgrid.dt;
                rf_prebf(:,:,aa) = sensor_data';
            end
            file_out = ['PWNE',num2str(round(tone_burst_freq/1e6)), ...
                'MHz_sam_att0p1inc0p',num2str(round(attInc)), ...
                'f',num2str(alpha_power*10),...
                '_BA6inc12_nc', num2str(tone_burst_cycles),...
                '_',num2str(source_strength/1000),'kPa'];
            disp(file_out)
            save(file_out,'rf_prebf','fs','c0','param'); clear rf_prebf;


            clear transducer
        end
    end
end

return