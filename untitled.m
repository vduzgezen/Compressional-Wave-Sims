clear all;

T_type = 1;
Ni = 1;
Di = 50;
Dhi = 500;
Dvi = 12;
Nl = 4;
Depth1 = 300;
Depth2 = 200; 
Depth3 = 300;
Depth4 = 1;
restrict = 1;
Rho_1 = 1000;
Rho_2 = 1300;
Rho_3 = 1600;
Rho_4 = 2000;
Rho_bt = 1000;
a_1 = 0.5;
a_2 = 0.7;
a_3 = 0.8;
a_4 = 3;
a_bt = 0.5;
BonA_1 = 5;
BonA_2 = 7;
BonA_3 = 8;
BonA_4 = 4;
BonA_bt = 5;
Fc = 4000000;
Coni = -20;
Foc_dep = 4000;
Bw = 4;
Ny_z = 108;
s1 = 1540;
s2 = 1600;
s3 = 1660;
s4 = 1700;

comp_wave_simulation_function(T_type, Nl, Ny_z, Depth1, Depth2, Depth3, Depth4, s1, s2, s3, s4,...
        Rho_1, Rho_2, Rho_3, Rho_4, Ni, Di, Dhi, Dvi, Coni, Foc_dep, Fc, Bw, restrict,...
        a_1, a_2, a_3, a_4, BonA_1, BonA_2, BonA_3, BonA_4)
    
    
    