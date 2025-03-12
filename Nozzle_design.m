clear; clc; close all;
format long

%Properties in the exit of the diffuser
T3   = 701.21121;           %Temperature     
                            %(T)   [K]
P3   = 189115.270458306;    %Pressure  
                            %(P)   [Pa]                         
M3   = 0.656150437410476;   %Mach Number                        
                            %(M)
rho3 = 0.939541230625329;   %Density                         
                            %(rho) [kg/m^3]
D    = 100/1000;            %External Diameter                  
                            %(D)   [mm]
Din  =(100-24.368)/1000;    %Internal Diameter                  
                            %(Din) [mm]
d3   =(D-Din);              %Diameter                           
                            %(d)   [m]
A3   = d3*1;                %Area for unit span wedge surface   
                            %(A)   [m^2]

g    = 1.41203;             %Degrees of freedom in a gas (gamma)            
R    = 287;                 %Gas Universal Constant [J*Kg^-1*K^-1]

a3   = sqrt(g*R*T3);        %Sound Speed [m/s]
u3   = a3*M3;               %Flow Velocity [m/s]

%Throat Area Auxiliar -> To evaluate the Mach Number after Expansion

Aaux   = A3/sqrt((1/M3^2)*(((2 /(g+1))*...
      (1+((g-1)/2)*M3^2))^((g+1)/(g-1))));
  
  
%Stagnation Properties in Volume 3

T0aux   = T3*(1+((g-1)/2)*M3^2);
P0aux  = P3/(1+((g-1)/2)*M3^2)^(-g/(g-1));
rho0aux = rho3*((1+((g-1)/2)*M3^2)^(1/(g-1)));

%Properties in Volume 4 (Total properties stay the same)



A4   = 0.1;

f    = @(M4) ((A4/Aaux)^2 - 1/M4^2 * (2/(g+1) * (1 + (g-1)...
       * M4^2 /2))^( (g+1)/(g-1))); 
M4   = fzero(f,0.01,3);
T4   = T0aux/(1+((g-1)/2)*M4^2);
P4   = P0aux*((1+(g-1)/2)*M4^2)*(-g/(g-1));
rho4 = rho0aux/((1+((g-1)/2)*M4^2)^(1/(g-1)));
a4   = sqrt(g*R*T4);
u4   = M4*a4;

%Heat addition (4-5)

P04  = P0aux;
T04  = T0aux;
T5   = T4+1000;

%Properties in Volume 5

A5 = A4;

f     = @(M5) T5/T4 - (((1+g*M4^2)/(1+g*M5^2))^2)*((M5/M4)^2);
M5    = fzero(f,M4);
P05   = P04* ((1+g*M4^2)/(1+g*M5^2))...
        *((1+((g-1)/2)*M5^2)/(1+((g-1)/2)*M4^2))^(g/(g-1));
T05   = T04*((1+g*M4^2)/(1+g*M5^2))...
        *((1+((g-1)/2)*M5^2)/(1+((g-1)/2)*M4^2))*(M5/M4)^2;
rho5  = rho4*(((1+g*M5^2)/(1+g*M4^2)))*(M4/M5)^2;
rho05 = rho5*(1+((g-1)/2)*M5^2)^(1/(g-1));
P5    = P4*((1+g*M4^2)/(1+g*M5^2));
a5    = sqrt(R*g*T5);
u5    = M5*a5;

%Properties in Volume 6

A6   = A5/sqrt((1/M5^2)*(((2/(g+1)))*(1+((g-1)/2)*M5^2))^((g+1)/(g-1)));

Dt = A6/1;                          %Area per unit span

%Calculating distances

x_t = sqrt((1.2*D+Dt)^2-(2*Dt+0.2*D)^2);
x_t2 = sqrt((1.5*D+Dt)^2-(2*Dt+0.5*D)^2)+x_t;


% Geometric Shapes

[xb,Rb_x,Mb_x,Pb_x,Tb_x,rhob_x,ub_x] = bell_geometry...
                                       (x_t,x_t2,Dt,P05,T05,rho05,A6);

[xc,Rc_x,Mc_x,Pc_x,Tc_x,rhoc_x,uc_x] = conical_geometry...
                                       (x_t,x_t2,Dt,P05,T05,rho05,A6);
                                   
% Comparisons

figure(13);
plot(xc,Rc_x,'k'); axis equal;
hold on
plot(xb,Rb_x,'r');
title('Geometry Comparison')
legend('Conical Nozzle','Bell-Shaped Nozzle');
ylim([0,0.11]);
xlim([0 xc(end)]);

figure(14)
plot(xc,Mc_x,'r');
hold on
plot(xb,Mb_x,'k');
title('Mach Comparison')
xlabel('Position (m)');
ylabel('Mach Number');
legend('Conical Nozzle','Bell-Shaped Nozzle');
xlim([0 xc(end)]);

figure(15);
plot(xc,Pc_x/Pc_x(1),'r');
hold on;
plot(xb,Pb_x/Pb_x(1),'k');
title('Pressure Ratio Comparison');
xlabel('Position (m)');
ylabel('Pressure ratio P/P0');
legend('Conical Nozzle','Bell-Shaped Nozzle');
xlim([0 xc(end)]);

figure(16);
plot(xc,Tc_x/Tc_x(1),'r');
hold on;
plot(xb,Tb_x/Tb_x(1),'k');
title('Temperature Ratio Comparison');
ylabel('Temperature ratio T/T0')
xlabel('Position(m)');
legend('Conical Nozzle','Bell-Shaped Nozzle');
xlim([0 xc(end)]);

figure(17);
plot(xc, rhoc_x,'r');
hold on
plot(xb, rhob_x,'k');
title('Density Comparison');
ylabel('Density [kg/m^3]');
xlabel('Position (m)');
legend('Conical Nozzle','Bell-Shaped Nozzle');
xlim([0 xc(end)]);

figure(18)
plot(xc,uc_x,'r');
hold on
plot(xb,ub_x,'k');
title('Flow Velocity Comparison');
xlabel('Position (m)');
ylabel('Velocity (m/s)');
legend('Conical Nozzle','Bell-Shaped Nozzle');
xlim([0 xc(end)]);