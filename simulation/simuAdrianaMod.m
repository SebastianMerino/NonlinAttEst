% clc
clear
clc
addpath(genpath(pwd))
addpath(genpath('/opt/MATLAB Add-Ons'))

alpha=[0.4 0.6 0.8];
ba = [6 9 12];
for rr=1:4
    for bb=1:length(ba)
        for al=1:length(alpha)
            % simulation settings
            DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
            RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

            % =========================================================================
            % DEFINE THE K-WAVE GRID
            % =========================================================================

            % set the size of the perfectly matched layer (PML)
            pml_x_size = 40;                % [grid points] *2
            pml_y_size = 10;                % [grid points]
            pml_z_size = 10;                % [grid points]

            Nx = 1248 - 2 * pml_x_size;      % [grid points]
            Ny = 675- 2 * pml_y_size;      % [grid points]
            Nz = 160 - 2 * pml_z_size;      % [grid points]

            % calculate the spacing between the grid points
            dx = 0.06e-3;                    % [m]
            x = dx*Nx;                      % [m]
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
            medium.alpha_coeff = alpha(al);      % [dB/(MHz^y cm)]
            medium.alpha_power = 1.05;

            % create the time array
            t_end = (Nx * dx) * 2.2 / c0;   % [s]
            kgrid.makeTime(c0, [], t_end);

            % =========================================================================
            % DEFINE THE INPUT SIGNAL
            % ========================================================================
            % define properties of the input signal
            tone_burst_freq = 5e6;        % [Hz]
            tone_burst_cycles = 3.5;
            source_strength_vector = [80 400]*1e3;

            for ss = 1:length(source_strength_vector)
                % append input signal used to drive the transducer
                source_strength =source_strength_vector(ss);          % [Pa]
                input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles, 'Envelope', 'Gaussian');
                input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

                % =========================================================================
                % DEFINE THE ULTRASOUND TRANSDUCER
                % =========================================================================

                % physical properties of the transducer
                transducer.number_elements = 128;  	% total number of transducer elements
                transducer.element_width = 5;       % width of each element [grid points] 9
                transducer.element_length = 100;  	% length of each element [grid points] 166
                transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points] 1
                transducer.radius = inf;            % radius of curvature of the transducer [m]

                % calculate the width of the transducer in grid points
                transducer_width = transducer.number_elements * transducer.element_width ...
                    + (transducer.number_elements - 1) * transducer.element_spacing;

                % use this to position the transducer in the middle of the computational grid
                transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

                % properties used to derive the beamforming delays
                transducer.sound_speed = c0;                    % sound speed [m/s]
                transducer.focus_distance = inf;              % focus distance [m]
                transducer.elevation_focus_distance = inf;    % focus distance in the elevation plane [m] 0.02
                transducer.steering_angle = 0;                  % steering angle [degrees]

                % apodization
                transducer.transmit_apodization = 'Rectangular';
                transducer.receive_apodization = 'Rectangular';

                % define the transducer elements that are currently active
                number_active_elements = 128;
                transducer.active_elements = ones(transducer.number_elements, 1);

                % append input signal used to drive the transducer
                transducer.input_signal = input_signal;

                % create the transducer using the defined settings
                transducer = kWaveTransducer(kgrid, transducer);

                % print out transducer properties
                transducer.properties;

                % =========================================================================
                % DEFINE THE MEDIUM PROPERTIES
                % =========================================================================

                % define a large image size to move across
                % number_scan_lines = 500;
                Nx_tot = Nx;
                Ny_tot=Ny;
                Nz_tot = Nz;
                BonA_map = ba(bb);
                alpha_coeff_map = alpha(al);

                medium.sound_speed_ref = 1540;
                medium.sound_speed = 1540;

                % =========================================================================
                % RUN THE SIMULATION
                % =========================================================================
                input_args = {...
                    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
                    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};


                matFile = ['BAPW_STD2_REF2024_',num2str(rr),'.mat'];
                medium.density = load(fullfile('random_media',matFile)).density_map;

                % update the command line status
                disp('');
                disp('Computing scan line PLANE WAVE...');

                medium.alpha_coeff = alpha_coeff_map;
                medium.BonA = BonA_map;


                % run the simulation
                sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});

                file_out = ['PWNE_5MHz_homo_BA',num2str(ba(bb)),...
                    '_att0p',num2str(alpha(al)*10,'%.1i'),'f1p2_nc3p5_',...
                    num2str(source_strength/1000),'kPa_sam',num2str(rr)];
                save(file_out,'sensor_data'); clear transducer,
            end
        end
    end
end