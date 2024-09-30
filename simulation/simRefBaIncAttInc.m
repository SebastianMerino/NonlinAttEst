% Simulating B-mode Ultrasound Images Example


clearvars; close all; clc;
%addpath(genpath([pwd,'/k-wave-toolbox-version-1.3']))
% simulation settings
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
%medium.alpha_coeff = 0.1;      % [dB/(MHz^y cm)]
medium.alpha_power = 2;
%medium.BonA = 11;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
% tone_burst_freq = 5e6;        % [Hz]

for tone_burst_freq_MHz = 4:6
    for ii = 1:3 % number of sims
        tone_burst_cycles = 10;

        source_strength_vector = [80 400]*1e3;
        for ss = 1:length(source_strength_vector)
            % append input signal used to drive the transducer
            source_strength =source_strength_vector(ss);          % [Pa]
            input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq_MHz*1e6, tone_burst_cycles);
            input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

            % =========================================================================
            % DEFINE THE ULTRASOUND TRANSDUCER
            % =========================================================================

            % physical properties of the transducer
            %transducer.number_elements = 32;  	% total number of transducer elements
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
            %number_active_elements = 32;
            transducer.active_elements = ones(transducer.number_elements, 1);



            transducer.input_signal = input_signal;
            % create the transducer using the defined settings
            transducer = kWaveTransducer(kgrid, transducer);


            % if ss == 1
            %     clear transducer
            %     continue
            % end
            % print out transducer properties
            transducer.properties;

            % =========================================================================
            % DEFINE THE MEDIUM PROPERTIES
            % =========================================================================

            % define a large image size to move across
            %number_scan_lines = 500;
            Nx_tot = Nx;
            Ny_tot = Ny;
            Nz_tot = Nz;

            %% Maps of B/A and AC
            switch ii
                case 1
                    att = 0.1; ba = 12;
                case 2
                    att = 0.1; ba = 6;
                case 3
                    att = 0.15; ba = 10;
            end
            att_map = single(att*ones(Nx,Ny_tot,Nz));
            BonA_map = single(ba*ones(Nx,Ny_tot,Nz));


            %%

            medium.sound_speed_ref = 1540;
            medium.sound_speed = 1540;
            % =========================================================================
            % RUN THE SIMULATION
            % =========================================================================

            % preallocate the storage
            %scan_lines = zeros(number_scan_lines, kgrid.Nt);

            % set the input settings
            input_args = {...
                'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
                'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

            % run the simulation if set to true, otherwise, load previous results from
            % disk
            %if RUN_SIMULATION

            % set medium position
            medium_position = 1;


            % loop through the scan lines
            for aa = 1:4% number slices
                % for scan_line_index = 1% 1:number_scan_lines
                filename = fullfile(pwd,'random_media',['BAPW_STD2_REF2024_',num2str(aa)]);
                %density = 1000;
                load(filename);
                density_map = density; clear density;
                % density_map = single(1000 * (1 + 0.02*randn(Nx,Ny_tot,Nz)));
                
                % update the command line status
                disp('');
                %disp(['Computing scan line ' num2str(scan_line_index) ' of ' num2str(number_scan_lines)]);
                disp(['Computing scan line PLANE WAVE...']);

                % load the current section of the medium
                %medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny - 1, :);
                medium.BonA = BonA_map;
                %medium.alpha_coeff = alpha_coeff_map(:, medium_position:medium_position + Ny - 1, :);
                medium.alpha_coeff = att_map;
                %medium.density = density_map(:, medium_position:medium_position + Ny - 1, :);
                medium.density = density_map;
                % run the simulation
                %sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});

                param.c0 = c0;
                param.att_map_setting = att_map(:,:,1);
                param.BonA_map_setting = BonA_map(:,:,1);
                param.att_power_law = medium.alpha_power;

                sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});

                %save FULLpwnoelev_sensor_data sensor_data;
                %save FULLpwnoelev_all_data;
                fs = 1/kgrid.dt;
                fnumber = 3;
                %rf(:,:,aa) = bf_planewave(sensor_data, fs, fnumber);
                rf_prebf(:,:,aa) = sensor_data';
            end

            file_out = ['PWNE',num2str(tone_burst_freq_MHz),'MHz_refBA',...
                num2str(ba),'_att0p',num2str(att*100,'%.2i'),...
                'f2_nc',num2str(tone_burst_cycles),'_',...
                num2str(source_strength/1000),'kPa'];
            save(file_out,'rf_prebf','fs','c0','param'); clear rf; clear transducer

        end




    end
end
return

