function oxy_ode = Oxychlorination_reactor_ode_team6(X,Y)

%%Output Variables
%Chemical Species Molar Flowrates
R = Y(1)                        ; % C2H4  (mol/L-hr)
T = Y(2)                        ; % HCl    (mol/L-hr)
U = Y (3)                       ; % O2     (mol/L-hr)
V = Y (4)                       ; % C2H4Cl2 (mol/L-hr)
W = Y (5)                       ; % H2O      (mol/L-hr)
Xx = Y (6)                      ; % C2H3Cl3  (mol/L-hr)
Yy = Y (7)                      ; % CO2      (mol/L-hr)
Z = Y (8)                       ; % Cl2      (mol/L-hr)


%Other Solution

Temp_cf = Y(9)                  ; % (K)  Temperature of the Coolant Fluid
Temp = Y(10)                    ; % (K)  Temperature in the Reactor
P = Y(11)                       ; % (kPa) (kPa) Pressure in the Reactor


%%Internal Variables
%Defing Constants
 
d_tube = 0.0254                 ; %(m)    Outer Tube Diameter (TA recommended this)                    
n_tubes = 294                   ; %       Number of Tubes
crossarea_tube = .149           ; %(m^2) Cross Sectional Area of Tube
bed_void = 0.50                 ; %      Voidage in the Bed
d_part = 0.003175               ; %(m)    Particle Diameter 

massflow_cf = 11.1111           ; %(kg/s)         Mass Flowrate of Coolant Fluid (Water)
MW_cf = 18.02                   ; %(Kg/kmol)      Molecular Weight of Coolant Fluid (Water)
Cp_cf = 4.18                    ; %(kj/mol-k)     Heat Capacity of Coolant Fluid (Water)
U_cf = 850                      ; % (W/m^2-k)     Heat Transfer Coeffiecient for Coolant Fluid (water)
Temp_out =533                   ; % (k)           Temperature Out
Gas_con = 0.008314              ; % (kJ/mol-k)    Ideal Gas Constant
U_a = 124.9                     ; % (W/m^2-k)     Heat transfer coefficient
Temp_a = 298                    ; %(K)


v_s = 0.0013944                 ; %(m/s)          Superlatice Gas Velocity 
void_frac = 0.5                 ; %                void fraction
gas_vis = 0.02                  ; % (cP)           Gas Viscocity
rho_int =  3950                 ; % (kg/m^3)       Particle Desntity

%Heats of Reaction 

Hr1 = -238.6148165              ; % (Kj/mol) Heat of Reaction for Reaction 1
Hr2 = -163.7605385              ; % (Kj/mol) Heat of Reaction for Reaction 2
Hr3 =-1322.753819               ; % (Kj/mol) Heat of Reaction for Reaction 3
Hr4f=-116.0417248               ; % (Kj/mol) Heat of Reaction for Reaction 4 (Forward)
Hr4r=116.0417248                ; % (Kj/mol) Heat of Reaction for Reaction 4 (Reverse)

% Feed Values
R_feed = 970.24                 ;%(Kmol/hr)
T_feed = 2910.75                ;%(Kmol/hr)
U_feed = 970.24                 ;%(Kmol/hr)
Z_feed = 5.82                   ;%(Kmol/hr)
Temp_cfin = 306                 ;% (K)
Temp_in = 298                   ;% (k)
P_in = 2000                     ;% (kPa)

%%ODE Definitions
%Species Heat capacities (@ 298K)

Cp_R = 32.083 + (-1.4831 * 10^(-2)) * Temp + (2.4774 * 10 ^ (-4)) * Temp^2    ; % Heat capacity of C2H4
Cp_T = 29.244 + (-1.2615 * 10^(-2)) * Temp + (1.1210 * 10 ^ (-6)) * Temp^2    ; % Heat capacity of HCl
Cp_U = 29.526 + (-8.8999 * 10^(-3)) * Temp + (3.8083 * 10 ^ (-5)) * Temp^2    ; % Heat capacity of O2
Cp_V = 15.730 + (-2.6124 * 10^(-1)) * Temp + (2.1489 * 10 ^ (-4)) * Temp^2    ; % Heat capacity of C2H4Cl2
Cp_W = 33.933 + (-8.4186 * 10^(-3)) * Temp + (2.9906  * 10 ^ (-5)) * Temp^2    ; % Heat capacity of H2O
Cp_Xx = 28.881 + (2.4893 * 10^(-1)) * Temp + (-3.4963 * 10 ^ (-4)) * Temp^2    ; % Heat capacity of C3H3Cl3
Cp_Yy = 27.437 + (4.2315 * 10^(-2)) * Temp + (-1.9555  * 10 ^ (-5)) * Temp^2    ; % Heat capacity of CO2
Cp_Z = 27.213 + (3.0426 * 10^(-2)) * Temp + (-3.335  * 10 ^ (-5)) * Temp^2    ; % Heat capacity of Cl2

%Arrhenius Reaction constants 

k1 = 10^4.2 * exp(-40.1/(0.008314 * Temp))        ;
k2 = 10^13.23 * exp(-128.04/(0.008314 * Temp))    ;
k3 = 10^6.78 * exp (17.13-13000/(1.987 * Temp))   ;
k4f = 1000*exp(17.13-1300/(1.987*Temp))           ;
k4r = exp(5.4 + 16000/(1.987*Temp))               ;
k4 = k4f/k4r                                      ;


%species partial Pressure (Assume Ideal gas law: p_i = (R * T * n_i)/ V)
p_R = Gas_con * Temp * R           ;
p_U = Gas_con * Temp * U           ;
p_V = Gas_con * Temp * V           ;
p_Z = Gas_con * Temp * Z           ;


%Reaction Rates

r1 = -k1 * p_R * (p_Z^ 0.5)            ;
r2 = -k2 * p_V * (p_Z ^ 0.5)           ;
r3 = -k3 * p_U  * (p_Z ^ 0.5) *p_R     ;
r4f = 0.5 * k4 * p_U * (p_Z ^ (-1))    ;
r4r = -0.5 * k4 * p_U * (p_Z ^ (-1))   ;

%Mole Balance Differentials
rR = R_feed - r1 -r3                           ;
rT = T_feed  - 2*r1 -r2 - 4*r4f + 4*r4r        ;
rU = U_feed - 0.5*r1 -0.5*r2 -3*r3 -r4f +r4r   ; 
rV = r1 - r2                                   ;
rXx = r2                                       ;
rYy = 2*r3                                     ;
rZ = Z_feed + 2*r4f - 2*r4r                    ;
rW =  r1 + r2 +2*r3 +2*r4f-2*r4r               ;


%Coolant Temp Differential 

rTemp_cf = (U_cf*(Temp_cf-Temp_a)) / ((massflow_cf / MW_cf) * Cp_cf) ;


% Thermal Differential 
sumHrxn = r1 * Hr1 + r2*Hr3 + r4f*Hr4f + r4r*Hr4r ;
sumCp = rR*Cp_R + rT*Cp_T + rU*Cp_U +rV*Cp_V + rW*Cp_W + rXx*Cp_Xx + rYy*Cp_Yy + rZ*Cp_Z ;
rTemp = (U_a*(Temp_a - Temp) + sumHrxn) / sumCp ;




%Pressure differential 

rP = (v_s/((d_part ^ 2) * (void_frac^3) * (crossarea_tube))) * (150 * gas_vis * (1-void_frac) + 1.5 * (v_s) * (rho_int) * (d_part)) ;

dR = rR           ;                                         
dT = rT           ;
dU = rU           ;
dV = rV           ;
dW = rW           ;
dXx = rXx         ;
dYy = rYy         ;
dZ = rZ           ;
dTemp_cf = rTemp_cf  ;
dTemp = rTemp       ; % (k)
dP = rP           ;

% Returens vertical ODE aray
oxy_ode = [dR; dT; dU; dV; dW; dXx; dYy; dZ; dTemp_cf; dTemp; dP] ;
end

