% clc

for al=2
     for cc=1
        addpath(genpath(pwd))
        addpath(genpath('/opt/MATLAB Add-Ons'))
        
        
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
        %alpha=[1 0.9 0.3 0.5 0.7 0.4 0.8 0.6];
	alpha=[0.08 0.12 0.16 0.2 0.18 0.06 0.1 0.14]
        cycles= [9 4.5 15 10];
        % define the properties of the propagation medium
        c0 = 1540;                      % [m/s]
        rho0 = 1000;                    % [kg/m^3]
        medium.alpha_coeff = alpha(al);      % [dB/(MHz^y cm)]
        medium.alpha_power = 2;
    %     medium.BonA = 6;
        % medium.alpha_mode='no_dispersion';
    
        % create the time array
        t_end = (Nx * dx) * 2.2 / c0;   % [s]
        kgrid.makeTime(c0, [], t_end);
    
        % =========================================================================
        % DEFINE THE INPUT SIGNAL
        % ========================================================================
        % define properties of the input signal
        % source_strength = 1e6;          % [Pa]
        tone_burst_freq = 5e6;        % [Hz]
        tone_burst_cycles = cycles(cc); %5 10 15 20  
        for vv=1
		v=[2 5 10 20]
       		 source_strength_vector = [40 10 20 50 100 200 80 400 800 500 1000 2000 4000 8000 5000 10000 20000]*1e3;
     
       		 for ss = 1:17
           	 % append input signal used to drive the transducer
           	 source_strength =source_strength_vector(ss);          % [Pa]
           	 input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles, 'Envelope', 'Gaussian');
           	 input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;
        
         	 % Window = tukeywin(length(input_signal),0.25); 
         	 % input_signal = input_signal.*Window';
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
        % Ny_tot = Ny + number_scan_lines * transducer.element_width;
        Ny_tot=Ny;
        Nz_tot = Nz;
    %     
         BonA_map = single(ones(Nx,Ny_tot,Nz));
             radius_disk = 9e-3;
             center_depth = 22.5e-3;
%     %         for mm=1:Nz
%     %             
%     %             BonA_map(:,:,mm) =single(randi([1 5],1) * makeDisc(Nx_tot, Ny_tot, round(center_depth/dx), Ny_tot/2, round(radius_disk/dx)));
%     %         end
%     %         BonA_map = BonA_map +randi([1 7],1);
       %      for mm=1:Nz
                 
      %           BonA_map(:,:,mm) =single(3* makeDisc(Nx_tot, Ny_tot, round(center_depth/dx), Ny_tot/2, round(radius_disk/dx)));
       %      end
    
        %    BonA_map = BonA_map + 6;
	 BonA_map = 6;

       % alpha = single(ones(Nx,Ny_tot,Nz));

	%alpha = single(ones(Nx,Ny_tot,Nz));
   	%	for mm=1:Nz
                 
         %        alpha(:,:,mm) =single(0.3* makeDisc(Nx_tot, Ny_tot, round(center_depth/dx), Ny_tot/2, round(radius_disk/dx)));
          %   end
    
           % alpha = alpha + 0.6;
	   %medium.alpha_coeff = alpha;

             alpha_coeff_map = alpha(al);
        %% 
        
        
        	medium.sound_speed_ref = 1540;
        	medium.sound_speed = 1540;
        
        
        
        % =========================================================================
        % RUN THE SIMULATION
        % =========================================================================
        	input_args = {...
            'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
            'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};
        
        % run the simulation if set to true, otherwise, load previous results from
        % disk
        % if RUN_SIMULATION
        
            % set medium position
            	medium_position = 1;
          	for aa = 33 % number slices
                	% for scan_line_index = 1% 1:number_scan_lines
                	filename = [pwd,'/density_map',num2str(aa)];
                	load(filename,'density')
                	density_map= density; clear density;
                % update the command line status
                	disp('');
                	disp(['Computing scan line PLANE WAVE...']);
    
                	medium.alpha_coeff = alpha_coeff_map;
                	medium.BonA = BonA_map;
                
                	medium.density = density_map;
                % run the simulation
                	sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});
                
            	end
		%file_out = ['RFinc_',num2str(source_strength/1000),'kPa_BA6-9_att06-09-33'];

           	file_out = ['RF9_',num2str(source_strength/1000),'kPa_BA6_att0',num2str(alpha(al)*10),'f2-33'];
	%	file_out = ['Bmode_hetero']
            	save(file_out,'sensor_data'); clear transducer
	    end
        end
    end
end

    
  
      

