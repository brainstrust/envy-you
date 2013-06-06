%% General constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Farad       = 96500         ;% [C mol-1]      Faradays constant
R_gas       = 8.315         ;% [J mol-1K-1]
Temp        = 300           ;% [K]

unitcon    = 10^(3)         ;% [-]            Factor to convert equations to another unit

%% Constants of Ostby, the values come from CELL ML

L_p         = 2.1e-9        ;% [m uM-1s-1]
R_tot       = 8.79e-8       ;% [m]
X_k         = 12.41e-3      ;% [uMm]
z_Na        = 1             ;% [-]
z_K         = 1             ;% [-]
z_Cl        = -1            ;% [-]
z_NBC       = -1            ;% [-]
g_K_k       = 40            ;% [ohm-1m-2]
g_KCC1_k    = 1e-2          ;% [ohm-1m-2]
g_NBC_k     = 7.57e-1       ;% [ohm-1m-2]
g_Cl_k      = 8.797e-1      ;% [ohm-1m-2]
g_NKCC1_k   = 5.54e-2       ;% [ohm-1m-2]
g_Na_k      = 1.314         ;% [ohm-1m-2]
J_NaK_max   = 1.42e-3       ;% [uMm s-1]
K_Na_k      = 10e3          ;% [uM]
K_K_s       = 1.5e3         ;% [uM]
k_C         = 7.35e-5       ;% [muM s-1]



%% Constants of the BK-channel

A_ef_k      = 3.7e-9        ;% m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
v_6         = 22e-3         ;% V        Constant
v_4         = 14.5e-3       ;% V        A measure of the spread of the distrubution
psi_w       = 2.664         ;% s-1      A characteristic time
G_BK_k      = 4.3e3         ; % pS      Constant estimation based on Ermentrout
g_BK_k      = G_BK_k*10^(-12)/A_ef_k ;% ohm-1m-2  Specific capacitance of the BK-Channel in units of Ostby
VR_pa       = 0.001       ;% [-]       The estimated volume ratio of perivascular space to astrocyte: Model estimation
VR_ps       = 0.001       ;% [-]       The estimated volume ratio of perivascular space to SMC: Model Estimation	




%% SMC constants


F_il = 7.5e2            ;%[-] scalingsfactor to fit the experimental data of Filosa
z_1 =4.5                ;%[-] parameter fitted on experimental data of Filosa
z_2 =-1.12e2            ;%[-] parameter fitted on experimental data of Filosa
z_3 =4.2e-1             ;%[-] parameter fitted on experimental data of Filosa
z_4 =-1.26e1            ;%[-] parameter fitted on experimental data of Filosa
z_5 =-7.4e-2            ;%[-] parameter fitted on experimental data of Filosa

% mVmicroM-1 The change in membrane potential by a scaling factor

% Koeningsberger et al.

Fmax_i		= 0.23;		% [microM/s]
Kr_i 		= 1; 		% [microM] ; Half saturation constant for agonist-dependent Ca$^{2+}$ entry
G_Ca		= 0.00129;	% [microM/mV/s]
v_Ca1		= 100;		% [mV]
v_Ca2		= -35; %-24;	% [mV]
R_Ca		= 8.5;		% [mV]
G_NaCa		= 0.00316;	% microM/mV/s
c_NaCa		= 0.5;		% microM
v_NaCa		= -30;
B_i			= 2.025;
cb_i		= 1;
C_i			= 55;
sc_i		= 2;
cc_i		= 0.9;
D_i			= 0.24;
vd_i		= -100;
Rd_i		= 250;
L_i			= 0.025;
gam			= 1970; % mVmicroM-1 The change in membrane potential by a scaling factor
F_NaK		= 0.0432;
G_Cl		= 0.00134;
v_Cl		= -25;
G_K			= 0.00446;
vK_i		= -94; 
lab 		= 45;
c_w			= 0;
bet			= 0.13;
v_Ca3		= -27;
R_K			= 12;
k_i			= 0.1;
K_d         = 1;            % = 1e3 nM Gonzalez
B_T         = 100;          % = 1e5 nM Gonzalez
Kinf_i      = 1e5;          % 100 mM K+ concentration in SMC

G_stretch   = 0.0061;       % uM mV-1 s-1
P_str       = 30;
Esac        = -18;          % mV
alpha1      = 0.0074;
sig0        = 500;


%% EC constants
% EC Koeningsberger et al.

Fmax_j		= 0.23;		% [microM/s]
Kr_j		= 1;
B_j 		= 0.5;
cb_j		= 1;
C_j			= 5;
sc_j		= 2;
cc_j		= 0.9;
D_j			= 0.24;
L_j			= 0.025;
G_cat 		= 0.66e-3;
E_Ca		= 50;
m3cat		= -0.18; %-6.18; %changed value!!! 
m4cat 		= 0.37;
J0_j 		= 0.029; %constant Ca influx (EC)
C_m 		= 25.8;
G_tot		= 6927;
vK_j 		= -80;
a1			= 53.3;
a2			= 53.3;
b			= -80.8;
c 			= -6.4; %-0.4; %changed value!!! 
m3b			= 1.32e-3;
m4b			= 0.3;
m3s			= -0.28;
m4s			= 0.389;
G_R			= 955;
v_rest		= -31.1;
k_j			= 0.1;





global CASE J_PLC

if CASE==1
g_hat 		= 50;
p_hat 		= 0;
p_hatIP3 	= 0.05;
elseif CASE==2
g_hat 		= 5;
p_hat 		= 0.05;
p_hatIP3 	= 0.05;
elseif CASE==3
g_hat 		= 0;
p_hat 		= 0;
p_hatIP3 	= 0.05;
elseif CASE==0
g_hat 		= 0;
p_hat 		= 0;
p_hatIP3 	= 0;
end

%% Myosin crossbridge model

K2_c        = 0.5;
K3_c        = 0.4;
K4_c        = 0.1;
K5_c        = 0.5;
K7_c        = 0.1;
gam_cross   = 17;

%% Koningsberger

% nu_r        = 100*inv(0.0075);
% sigp0_r     = 0.0191*inv(0.0075);
% kp_r        = 0.15;
% r0_r        = 20;
% siga0_r     = 1.8e5;
% ka_r        = 0.0006;
% ra_r        = 12;
% rb_r        = 15;
% hb_r        = 3;
% P_r         = 30*inv(0.0075);


%% Kelvin Voigt

P_r         = 4000;
rb_r        = 20e-6;
h0_r        = 3e-6;
R0pas_r     = 20e-6;
R0act_r     = 12e-6;
Etot_r      = 233e3;
Eact_r      = 233e3;
Epas_r      = 66e3;
nu_r        = 1e4;
