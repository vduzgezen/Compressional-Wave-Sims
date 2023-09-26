% ************************************************************************
% Description:
%   This script will take an excel file as input where the details of 
%   the excel file have been described in the template. This script runs as
%   many simulations as rows present in the excel file. Use the template
%   called Input_Template_Comp_Wave_Simulation.xlsx for accurate results.
% 
%   Vedat Can Duzgezen | MGH-CURT
% ************************************************************************

% Specify the file path of the Excel file (May vary)
file_path = 'Input_Template_Comp_Wave_Simulation.xlsx';

% Read the Excel file
data = xlsread(file_path);

% Get the number of rows and columns in the Excel file
[num_rows, num_cols] = size(data);

% Set the total number of simulations
total_simulations = num_rows - 2;

% Loop through each row and run the simulation
for i = 3:num_rows
    
    % Assign the columns of the current row to variables
    
    % Assign the The Tissue Type (Circular Lesions, Elliptical Lesions,
    % Layered Tissue, Layered + Circular, Layered + Elliptical)
    T_type = data(i, 1);
    if isnan(T_type) && isnan(data(2, 1))
        T_type = data(1, 1);
    elseif isnan(T_type)
        T_type = randi([data(1, 1), data(2, 1)], 1);
    end
    dataf(i, 1) = T_type;
    
    % Assign the Number of Layers
    Nl = data(i, 2);
    if isnan(Nl) && isnan(data(2, 2))
        Nl = data(1, 2);
    elseif isnan(Nl)
        Nl = randi([data(1, 2), data(2, 2)], 1);
    end
    dataf(i, 2) = Nl;
    
    % Assign the Phantom Depth
    Ny_z = data(i, 3);
    if isnan(Ny_z) && isnan(data(2, 3))
        Ny_z = data(1, 3);
    elseif isnan(Ny_z)
        Ny_z = randi([data(1, 3), data(2, 3)], 1);
    end
    dataf(i, 3) = Ny_z;
    
    % Assign the Thickness of Layer 1
    Depth1 = data(i, 4);
    if isnan(Depth1) && isnan(data(2, 4))
        Depth1 = data(1, 4);
    elseif isnan(Depth1)
        Depth1 = randi([data(1, 4), data(2, 4)], 1);
    end
    dataf(i, 4) = Depth1;
    
    % Assign the Thickness of Layer 2
    Depth2 = data(i, 5);
    if isnan(Depth2) && isnan(data(2, 5))
        Depth2 = data(1, 5);
    elseif isnan(Depth2)
        Depth2 = randi([data(1, 5), data(2, 5)], 1);
    end
    dataf(i, 5) = Depth2;

    % Assign the Thickness of Layer 3
    Depth3 = data(i, 6);
    if isnan(Depth3) && isnan(data(2, 6))
        Depth3 = data(1, 6);
    elseif isnan(Depth3)
        Depth3 = randi([data(1, 6), data(2, 6)], 1);
    end
    dataf(i, 6) = Depth3;
    
    % The thickness of layer 4 will be calculated from the left over space
    % in the  phantom, if the tissue is assigned to be 4 layered
    if T_type == 3 || T_type == 4 || T_type == 5 && Rho_4 == 3
        Depth4 = Ny_z - Coni - Dvi - Depth3;
    else
        Depth4 = NaN;
    end
    dataf(i, 7) = Depth4;
    
    % Assign the Speed of Sound for Layer 1
    s1 = data(i, 7);
    if isnan(s1) && isnan(data(2, 7))
        s1 = data(1, 7);
    elseif isnan(s1)
        s1 = randi([data(1, 7), data(2, 7)], 1);
    end
    dataf(i, 8) = s1;
    
    % Assign the Speed of Sound for Layer 2
    s2 = data(i, 8);
    if isnan(s2) && isnan(data(2, 8))
        s2 = data(1, 8);
    elseif isnan(s2)
        s2 = randi([data(1, 8), data(2, 8)], 1);
    end
    dataf(i, 9) = s2;
    
    % Assign the Speed of Sound for Layer 3
    s3 = data(i, 9);
    if isnan(s3) && isnan(data(2, 9))
        s3 = data(1, 9);
    elseif isnan(s3)
        s3 = randi([data(1, 9), data(2, 9)], 1);
    end
    dataf(i, 10) = s3;
    
    % Assign the Speed of Sound for Layer 4
    s4 = data(i, 10);
    if isnan(s4) && isnan(data(2, 10))
        s4 = data(1, 10);
    elseif isnan(s4)
        s4 = randi([data(1, 10), data(2, 10)], 1);
    end
    dataf(i, 11) = s4;
    
    % Assign the Density of Layer 1
    Rho_1 = data(i, 11);
    if isnan(Rho_1) && isnan(data(2, 11))
        Rho_1 = data(1, 11);
    elseif isnan(Rho_1)
        Rho_1 = randi([data(1, 11), data(2, 11)], 1);
    end
    dataf(i, 12) = Rho_1;
    
   % Assign the Density of Layer 2
    Rho_2 = data(i, 12);
    if isnan(Rho_2) && isnan(data(2, 12))
        Rho_2 = data(1, 12);
    elseif isnan(Rho_2)
        Rho_2 = randi([data(1, 12), data(2, 12)], 1);
    end
    dataf(i, 13) = Rho_2;

    % Assign the Density of Layer 3
    Rho_3 = data(i, 13);
    if isnan(Rho_3) && isnan(data(2, 13))
        Rho_3 = data(1, 13);
    elseif isnan(Rho_3)
        Rho_3 = randi([data(1, 13), data(2, 13)], 1);
    end
    dataf(i, 14) = Rho_3;
    
    % Assign the Density of Layer 4
    Rho_4 = data(i, 14);
    if isnan(Rho_4) && isnan(data(2, 14))
        Rho_4 = data(1, 14);
    elseif isnan(Rho_4)
        Rho_4 = randi([data(1, 14), data(2, 14)], 1);
    end
    dataf(i, 15) = Rho_4;
    
    % Assign the Number of Inclusions
    Ni = data(i, 15);
    if isnan(Ni) && isnan(data(2, 15))
        Ni = data(1, 15);
    elseif isnan(Ni)
        Ni = randi([data(1, 15), data(2, 15)], 1);
    end
    dataf(i, 16) = Ni;
    
    % Assign the Circular Inclusion Diameter
    Di = data(i, 16);
    if isnan(Di) && isnan(data(2, 16))
        Di = data(1, 16);
    elseif isnan(Di)
        Di = randi([data(1, 16), data(2, 16)], 1);
    end
    dataf(i, 17) = Di;
    
    % Assign the Elliptical Horizontal Diameter
    Dhi = data(i, 17);
    if isnan(Dhi) && isnan(data(2, 17))
        Dhi = data(1, 17);
    elseif isnan(Dhi)
        Dhi = randi([data(1, 17), data(2, 17)], 1);
    end
    dataf(i, 18) = Dhi;
    
    % Assign the Elliptical Vertical Diameter
    Dvi = data(i, 18);
    if isnan(Dvi) && isnan(data(2, 18))
        Dvi = data(1, 18);
    elseif isnan(Dvi)
        Dvi = randi([data(1, 18), data(2, 18)], 1);
    end
    dataf(i, 19) = Dvi;
    
    % Assign the Contrast of Inclusions
    Coni = data(i, 19);
    if isnan(Coni) && isnan(data(2, 19))
        Coni = data(1, 19);
    elseif isnan(Coni)
        Coni = randi([data(1, 19), data(2, 19)], 1);
    end
    dataf(i, 20) = Coni;
    
    % Assign the Focal Depth
    Foc_dep = data(i, 20);
    if isnan(Foc_dep) && isnan(data(2, 20))
        Foc_dep = data(1, 20);
    elseif isnan(Foc_dep)
        Foc_dep = randi([data(1, 20), data(2, 20)], 1);
    end
    dataf(i, 21) = Foc_dep;
    
    % Assign the Center Frequency
    Fc = data(i, 21);
    if isnan(Fc) && isnan(data(2, 21))
        Fc = data(1, 21);
    elseif isnan(Fc)
        Fc = randi([data(1, 21), data(2, 21)], 1);
    end
    dataf(i, 22) = Fc;
    
    % Assign Signal Bandwith
    Bw = data(i, 22);
    if isnan(Bw) && isnan(data(2, 22))
        Bw = data(1, 22);
    elseif isnan(Bw)
        Bw = randi([data(1, 22), data(2, 22)], 1);
    end
    dataf(i, 23) = Bw;
    
    % Assign the restriction coefficient which, if not NaN, in a multi
    % layered tissue, restricts the circular or elliptical lesions to a
    % singe layer ( Ex: restrict = 1 would make so all lesions would be
    % contained in the top layer)
    restrict = data(i, 23);
    if isnan(restrict) && isnan(data(2, 23))
        restrict = data(1, 23);
    elseif isnan(restrict)
        restrict = randi([data(1, 23), data(2, 23)], 1);
    end
    dataf(i, 24) = restrict;
    
    % Assign the Speed of Sound in Soft Inclusions
    a_1 = data(i, 24);
    if isnan(a_1) && isnan(data(2, 24))
        a_1 = data(1, 24);
    elseif isnan(a_1)
        a_1 = randi([data(1, 24), data(2, 24)], 1);
    end
    dataf(i, 25) = a_1;

    % Assign the Speed of Sound in Hard Inclusions
    a_2 = data(i, 25);
    if isnan(a_2) && isnan(data(2, 25))
        a_2 = data(1, 25);
    elseif isnan(a_2)
        a_2 = randi([data(1, 25), data(2, 25)], 1);
    end
    dataf(i, 26) = a_2;
    
    % Assign the Attenuation Coefficient of Layer 1
    a_3 = data(i, 26);
    if isnan(a_3) && isnan(data(2, 26))
        a_3 = data(1, 26);
    elseif isnan(a_3)
        a_3 = data(1, 26) + rand * (data(2, 26) - data(1, 26));
    end
    dataf(i, 27) = a_3;
    
    % Assign the Attenuation Coefficient of Layer 2
    a_4 = data(i, 27);
    if isnan(a_4) && isnan(data(2, 27))
        a_4 = data(1, 27);
    elseif isnan(a_4)
        a_4 = data(1, 27) + rand * (data(2, 27) - data(1, 27));
    end
    dataf(i, 28) = a_4;
    
    % Assign the Attenuation Coefficient of Layer 3
    BonA_1 = data(i, 28);
    if isnan(BonA_1) && isnan(data(2, 28))
        BonA_1 = data(1, 28);
    elseif isnan(BonA_1)
        BonA_1 = data(1, 28) + rand * (data(2, 28) - data(1, 28));
    end
    dataf(i, 29) = BonA_1;
    
    % Assign the Attenuation Coefficient of Layer 4
    BonA_2 = data(i, 29);
    if isnan(BonA_2) && isnan(data(2, 29))
        BonA_2 = data(1, 29);
    elseif isnan(BonA_2)
        BonA_2 = data(1, 29) + rand * (data(2, 29) - data(1, 29));
    end
    dataf(i, 30) = BonA_2;
    
    % Assign the Nonlinearity Parameter of Layer 1
    BonA_3 = data(i, 30);
    if isnan(BonA_3) && isnan(data(2, 30))
        BonA_3 = data(1, 30);
    elseif isnan(BonA_3)
        BonA_3 = data(1, 30) + rand * (data(2, 30) - data(1, 30));
    end
    dataf(i, 31) = BonA_3;
    
    % Assign the Nonlinearity Parameter of Layer 2
    BonA_4 = data(i, 31);
    if isnan(BonA_4) && isnan(data(2, 31))
        BonA_4 = data(1, 31);
    elseif isnan(BonA_4)
        BonA_4 = data(1, 31) + rand * (data(2, 31) - data(1, 31));
    end
    dataf(i, 32) = BonA_4;
    
    % Display the progress
    fprintf('\n==> Running Simulation %d out of %d... \n', i - 2 , total_simulations);
    
    % Run the simulation function with the assigned variables
    comp_wave_simulation_function(T_type, Nl, Ny_z, Depth1, Depth2, Depth3, Depth4, s1, s2, s3, s4,...
        Rho_1, Rho_2, Rho_3, Rho_4, Ni, Di, Dhi, Dvi, Coni, Foc_dep, Fc, Bw, restrict,...
        a_1, a_2, a_3, a_4, BonA_1, BonA_2, BonA_3, BonA_4)
end

% Upload the determined inputs to an excel file for future reference after
% naming the table that contains the values.
inputFile = 'determined_inputs_comp_wave_simulation.xlsx';
variableNames = {'Tissue Type Code (Refer to Table)',...
    'Number of Layers', 'Phantom depth (mm)',...
    'Layer 1 Depth (mm)', 'Layer 2 Depth (mm)', 'Layer 3 Depth (mm)', 'Layer 4 Depth (mm)',...
    'Layer 1 Speed of Sound (m/s)', 'Layer 2 Speed of Sound (m/s)',...
    'Layer 3 Speed of Sound (m/s)', 'Layer 4 Speed of Sound (m/s)',...
    'Density of Layer 1', 'Density of Layer 2', 'Density of Layer 3',...
    'Density of Layer 4', 'Number of Inclusions', 'Diameter of Circular Inclusions (mm)',...
    'Horizontal Diameter of Elliptical Inclusions (mm)', 'Vertical Diameter of Elliptical Inclusions (mm)',...
    'Contrast of Inclusions (dB)', 'Focal Depth (mm)', 'Center Frequency (MHz)', 'Signal Bandwidth (%)',...
    'Inclusions Restricted into Layer #', 'Speed of Sound in Soft Inclusions (m/s)', 'Speed of Sound in Hard Inclusions (m/s)'};
dataTable = array2table(dataf(3:end,1:32), 'VariableNames', variableNames);
writetable(dataTable, inputFile, 'Sheet', 1, 'Range', 'A1');

