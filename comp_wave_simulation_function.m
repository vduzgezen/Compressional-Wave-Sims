function [] = comp_wave_simulation_function(T_type, Nl, Ny_z, Depth1, Depth2, Depth3, Depth4, s1, s2, s3, s4,...
        Rho_1, Rho_2, Rho_3, Rho_4, Ni, Di, Dhi, Dvi, Coni, Foc_dep, Fc, Bw, restrict,...
        a_1, a_2, a_3, a_4, BonA_1, BonA_2, BonA_3, BonA_4)
% Simulating ultrasound pressure field for a point source transmit
% 2-D
%
% based on the example for a B-mode scan using a linear array
%
% last updated on 01/11/2023
%
% Vedat Duzgezen & Marko Jakovljevic | MGH Curt

% simulation settings
RUN_SIMULATION = false;
VISUAL = false;
DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 40;                % [grid points]
pml_y_size = 40;                % [grid points]

% set desired grid size in 2-D not including the PML
x = Ny_z * 1e-3;                 % [m]
y = 128e-3;                      % [m]

% set the spacing between the grid points
dx = 1e-3;                       % [m]
dy = dx;                         % [m]

% set total number of grid points not including the PML
Nx = round(x/dx);
Ny = round(y/dy);

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 2500;                      % [m/s]
rho0 = 1900;  
% c0 = 1540;                      % [m/s]
% rho0 = 1000;                    % [kg/m^3]

% create the time array
% CFL = 0.2; % stability criteria for convergence
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
% tone_burst_freq = params.f0;        % [Hz]
tone_burst_freq = Fc; % Center Frequency
tone_burst_cycles = 1 / Bw; % 1/Band Width

P.tone_burst_freq = tone_burst_freq;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 64;  	% number of elements per subapeture: HARD-CODED
pitch = 0.3e-3;
transducer.element_width = round(pitch/dy);       % width of each element [grid points]
transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

transducer.transducer_width = transducer_width;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, round(Ny/2 - transducer_width/2)]);

transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = Foc_dep * 1e-3;              % focus distance [m]

transducer.input_signal = input_signal;

N_points = transducer.number_elements * transducer.element_width; % not accounting for kerf

% Rx aperture
sensor.mask = zeros(Nx,Ny);
idc_y = (1:N_points) + transducer.position(2);

x_axis = (0:Nx-1) * dx;
Rx_depth = 67e-3;
% Rx_depth = 64e-3;
idc_x = find(x_axis>=Rx_depth,1,'first');

sensor.mask(idc_x,idc_y) = 1;

% sensor.record = {'p_final', 'p_max', 'p_rms'};
sensor.record = {'p_max_all','p'};

% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

Nx_tot = Nx;
Ny_tot = Ny;

P.Ny_tot = Ny_tot;
P.Nx_tot = Nx_tot;

%%

% Background Specle Parameters
alpha_power_map = 1.5;
background_map_mean = 1;
background_map_std = 0.004;
inclusion_indices = zeros(Nx, Ny);
%%

% Define the mesh grid for lesion formation
[nnx,nny]=meshgrid(1:Ny, 1:Nx);

% Define the shear sound speed [m/s]
sound_speed_map = ones(Nx, Ny) * s1;

if T_type == 3 || T_type == 4 || T_type == 5

    % LAYERED TISSUE MODEL

    % define the shear sound speed [m/s]
    sound_speed_map = ones(Nx, Ny) * s1;

    % Assign elasticities to respective layers
    if Nl == 2
        sound_speed_map(Depth1 + 1 : end, :) = s2;
        Depth3 = 0;
        Depth4 = 0;
    end
    if Nl == 3
        % Calculate the thickness of the third layer
        sound_speed_map(Depth1 + 1 : Depth1 + Depth2 , :) = s2;
        Depth3 = Nx - Depth1 - Depth2;
        sound_speed_map(Depth1 + Depth2 + 1 : end , :) = s3;
        Depth4 = 0;
    end
    if Nl == 4
        % Calculate the thickness of the fourth layer
        sound_speed_map(Depth1 + 1 : Depth1 + Depth2 , :) = s2;
        sound_speed_map(Depth1 + Depth2 + 1 : Depth1 + Depth2 + Depth3 , :) = s3;
        sound_speed_map(Depth1 + Depth2 + Depth3 + 1 : end , :) = s4;
        Depth4 = Nx - Depth1 - Depth2 - Depth3;
    end
else

    % If the tissue is not layered set the restriction to NaN
    Depth3 = 0;
    Depth4 = 0;
    restrict = NaN;

end

% Define mean & std for speckle targets
background_map = background_map_mean + background_map_std * randn([Nx_tot, Ny_tot]);

% Apply Speckles to the Layered Medium
sound_speed_map2 = sound_speed_map;
sound_speed_map = sound_speed_map .* background_map;

% Speckle Coefficient 

if T_type == 1 || T_type == 4

    % CIRCULAR

    % Set horizontal and vertical diameters equal to each other since the
    % elliptical incluson will be a circle under the conditions
    Dhi = Di;
    Dvi = Di;
end

if T_type == 1 || T_type == 2 || T_type == 4 || T_type == 5

    % ELLIPTICAL

    % Generate empty grid to track locations of soft inclusions
    soft_inclusions_x = zeros(Ni, 1);
    soft_inclusions_y = zeros(Ni, 1);

    % Set maximum duration for the loop
    max_duration = 10;  % Maximum duration in seconds
    start_time = tic();  % Start the timer

    for i = 1 : Ni
        % Generate random coordinates
        x = randi(Ny + round(Dhi/2));
        y = randi(Nx + round(Dvi/2));

        % Check for interference with existing inclusions and restrict to the specified layer if not NaN
        while any(sqrt(((soft_inclusions_x(1:i-1) - x)./Dhi).^2 + ((soft_inclusions_y(1:i-1) - y)./Dvi).^2) < 1) || ...
              any(sqrt(((soft_inclusions_x(1:i-1) - x)./Dhi).^2 + ((soft_inclusions_y(1:i-1) - y)./Dvi).^2) < 2) || ...
              (~isnan(restrict) && ((restrict == 1 && y > Depth1 - Dvi) || (restrict == 2 && (y < Depth1 + Dvi || y >= Nx - Depth3 - Dvi)) ...
              || (restrict == 3 && (y < Depth1 + Depth2 + Dvi || y >= Nx - Depth4 - Dvi)) || (restrict == 4 && y < Depth1 + Depth2 + Depth3 + Dvi)))
            x = randi(Ny + round(Dhi/2));
            y = randi(Nx + round(Dvi/2));

            % Check if the maximum duration has been reached
            elapsed_time = toc(start_time);
            if elapsed_time >= max_duration
                break;  % Exit the loop if the maximum duration has been reached
            end
        end

        % Set inclusion coordinates
        soft_inclusions_x(i) = x;
        soft_inclusions_y(i) = y;
        
        % Calculate inclusion indices
        inclusion_indices2 = sqrt(((nnx - x)./Dhi).^2 + ((nny - y)./Dvi).^2) < 1;
        inclusion_indices = inclusion_indices + inclusion_indices2;
        
    end
    
        % Apply Inclusions
        inclusion_values = sound_speed_map2(inclusion_indices > 0);
        inclusion_std = background_map_std * 10^(Coni/20);
        sound_speed_map(inclusion_indices > 0) = inclusion_values .* (background_map_mean + inclusion_std * randn(size(inclusion_values > 0)));

end

if T_type == 3 || T_type == 4 || T_type == 5

    % LAYERED MEDIA DENSITY

    % Define the layered medium density [m/s]
    density_map = ones(Nx, Ny) * Rho_1;
    alpha_coeff_map = ones(Nx, Ny) * a_1;
    BonA_map = ones(Nx, Ny) * BonA_1;

    % Assign elasticities to respective layers
    if Nl == 2
        density_map(Depth1 + 1 : end, :) = Rho_2;
        alpha_coeff_map(Depth1 + 1 : end, :) = a_2;
        BonA_map(Depth1 + 1 : end, :) = BonA_2;
        
        
    end
    if Nl == 3
        % Assign the parameters of the three layers
        density_map(Depth1 + 1 : Depth1 + Depth2 , :) = Rho_2;
        density_map(Depth1 + Depth2 + 1 : end , :) = Rho_3;
        alpha_coeff_map(Depth1 + 1 : Depth1 + Depth2 , :) = a_2;
        alpha_coeff_map(Depth1 + Depth2 + 1 : end , :) = a_3;
        BonA_map(Depth1 + 1 : Depth1 + Depth2 , :) = BonA_2;
        BonA_map(Depth1 + Depth2 + 1 : end , :) = BonA_3;
    end
    if Nl == 4
        % Assign the parameters of the four layers
        density_map(Depth1 + 1 : Depth1 + Depth2 , :) = Rho_2;
        density_map(Depth1 + Depth2 + 1 : end , :) = Rho_3;
        density_map(Depth1 + Depth2 + Depth3 + 1 : end , :) = Rho_4;
        alpha_coeff_map(Depth1 + 1 : Depth1 + Depth2 , :) = a_2;
        alpha_coeff_map(Depth1 + Depth2 + 1 : end , :) = a_3;
        alpha_coeff_map(Depth1 + Depth2 + Depth3 + 1 : end , :) = a_4;
        BonA_map(Depth1 + 1 : Depth1 + Depth2 , :) = BonA_2;
        BonA_map(Depth1 + Depth2 + 1 : end , :) = BonA_3;
        BonA_map(Depth1 + Depth2 + Depth3 + 1 : end , :) = BonA_4;
    end
else
    % Define the mass default density [kg/m^3]
    density_map = a_1 * ones(Nx, Ny);
    alpha_coeff_map = a_1 * ones(Nx, Ny);
    BonA_map = BonA_1 * ones(Nx, Ny);
end

figure;
imagesc(sound_speed_map);
colorbar;
% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

dir_name = 'elipse_2D_FSA';

% run the simulation if set to true
if RUN_SIMULATION

    % load the current section of the medium        
    medium.sound_speed = sound_speed_map;
    medium.density = density_map;
    medium.alpha_coeff = alpha_coeff_map;      % [dB/(MHz^y cm)]
    medium.alpha_power = alpha_power_map;
    medium.BonA = BonA_map;

    % =============================
    % cycle through sources
    % =============================
    ele_id = round(N_points/2);

    % update the command line status
    disp(['Simulating data for source # ' num2str(ele_id) ' of ' num2str(N_points) ' points.']);

    % apod_idc = idc_y([1:transducer.element_width] + (ele_id-1)*transducer.element_width);
    apod_idc = round(mean(idc_y)); % a single pixel active source

    % set source location
    x_axis = (0:Nx_tot-1) * dx; Tx_depth = 20e-3;
    idc_x = find(x_axis>=Tx_depth,1,'first');
    
    source.p_mask = zeros(Nx,Ny); % reset source mask
    source.p_mask(idc_x,apod_idc) = 1; % apodize
    source.p = [repmat(input_signal,[length(apod_idc),1]) zeros(length(apod_idc),length(input_signal))]; % reset waveforms

    % plot Tx and Rx config
    if VISUAL
        figure;
        imagesc(y_axis*1e3, x_axis*1e3, sound_speed_map); colorbar
        xlabel('lateral (mm)')
        ylabel('depth (mm)')
        title('Sound speed (m/s)')
        axis image; grid on
        hold on
        plot(y_axis(apod_idc)*1e3,Tx_depth*1e3,'x r')
        plot(y_axis(idc_y)*1e3,Rx_depth*1e3,'x k')
        legend('Active sources', 'Receiver locations')
    end

    Run the simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
end

% =========================================================================
% Beamforming
% =========================================================================
% 
% p_max = sensor_data.p_max_all;
% p_max = p_max/max(p_max(:));
% kgrid.makeTime
% signal = x, y or sensor_data.p?
% focal points vector? 
% element psotion = transducer.position?
end

