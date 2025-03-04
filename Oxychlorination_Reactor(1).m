
clear
clc
clear figure

%%Volume Domain
v_tot = 6                          ; %(m^3)
v_domain = linspace(0,1,100)       ;

%%Catalyst Information
rho_bulk = 1975                    ; %(kg*cat/m^3)  Catalyt Bulck Density
bed_void = 0.50                    ; %             Voidage in the Bed
rho_int = 3950                     ; %(kg*cat/m^3) Intrinsic Catalyst Density


%%Constants
MW_cf = 18.02                      ; %(kg/kgmol) MW of Cooling Fluid
Tau = 3.6                          ; %(s)   Space Time
Temp_out = 533                     ; %(K) Temperature Out
Gas_con = 0.008314                 ; %(kj/mol-k) Ideal Gas Constant

%% Initial Conditions
R_in = 970.24                      ; %(kgmol/hr) C2H4 Feed Flowrate
T_in = 2910.75                     ; %(kgmol/hr) HCl
U_in = 970.24                      ; %(kgmol/hr) O2
V_in = 0                           ; %(kgmol/hr) No C2H4Cl2 in Feed Flowrate
W_in = 0                           ; %(kgmol/hr) No H2O in Feed
Xx_in = 0                          ; %(kgmol/hr) No C2H3Cl3 in Feed
Yy_in = 0                          ; %(kgmol/hr) No CO2 in Feed
Z_in = 5.82                        ; %(kgmol/hr) Cl2 feed Flowrate
Temp_cfin = 306                     ; %(K) Temp of cooling Fluid from Table 12.1
Temp_in = 298                      ; %(K) Temperature In
P_in = 2000                        ; %(kPa) Inlet Pressure

%%Collect Initial Conditions into an Array

IC = [R_in T_in U_in V_in W_in Xx_in Yy_in Z_in Temp_cfin Temp_in P_in] ;

%% Solve ODE
[Xsol, Ysol] = ode45(@Oxychlorination_reactor_ode_team6, v_domain, IC);


%%Data Handling -Extracting Solution for Each Variable

IC = [R_in T_in U_in V_in W_in Xx_in Yy_in Z_in Temp_cfin Temp_in P_in] ;
R_sol = Ysol(:,1) ;
T_sol = Ysol(:,2) ;
U_sol = Ysol(:,3) ;
V_sol = Ysol(:,4) ; 
W_sol = Ysol (:,5);
Xx_sol = Ysol (:,6);
Yy_sol = Ysol(:,7);
Z_sol = Ysol (:,8) ;
Temp_cf_sol = Ysol(:,9);
Temp_sol =Ysol(:,10);
P_sol = Ysol (:,11) ;


%%Plots

fig1 = figure(1)
%set(fig1, 'Carlos' , 'C2H4') ;
area(Xsol,R_sol)                 ;
xlabel('Volume (m^3)')           ;
ylabel('C2H4 Flowrate (kgml/hr)');

fig2 = figure(2)
%set(fig2, 'Carlos' ,'HCl') ;
area(Xsol,T_sol);
xlabel('Volume (m^3)')  ;
ylabel('HCl Flowrate (kmol/hr)');

fig3 = figure(3)
%set(fig3, 'Carlos' , '02')
area(Xsol, U_sol);
xlabel('Volume (m^3)');
ylabel('02 Flowrate (kmol/hr)');

fig4 = figure(4)
%set(fig3, 'Carlos', 'Reactor Temperature'
area(Xsol, Temp_sol);
xlabel('Volume (m^3)');
ylabel('Temperature (k)');

fig5 = figure(5)
%set(fig3, 'Carlos','Coolant Temperature')
area(Xsol, Temp_cf_sol);
xlabel('Volume (m^3)')
ylabel('Coolant Temp(k)');

fig6 = figure(6) 
%set(fig3, 'Carlos','Pressure')
area(Xsol,P_sol);
xlabel('Volume (m^3)');
ylabel('Pressure (kPa)');




