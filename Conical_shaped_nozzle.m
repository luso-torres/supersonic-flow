function [xc,Rc_x,Mc_x,Pc_x,Tc_x,rhoc_x,uc_x] = conical_geometry...
                                     (x_t,x_t2,Dt,P05,T05,rho05,A6)

g    = 1.41203;             %Degrees of freedom in a gas (gamma)            
R    = 287;                 %Gas Universal Constant [J*Kg^-1*K^-1]
D    = 0.1;                 %External Diameter

%Auxiliar function to round numbers
roundn = @(x,n) round(x.*10.^n)./10.^n; 

xc1 = 0:0.00001:x_t;
xc2 = x_t:0.00001:roundn(x_t2,5);
xc  = [xc1 xc2];

Rc_x1 = (Dt-D)/x_t*xc1+D;
Rc_x2 = (D-Dt)/(xc(end)-xc2(1))*(xc2 - xc2(1))+ Dt;
Rc_x  = [Rc_x1 Rc_x2];

Ac_x = Rc_x;
Mc_x = zeros(length(Ac_x),1);

%Subsonic Region

for i = 1:find(Ac_x == min(Ac_x))
    A_x = Ac_x(i);
    f = @(Mc_x) ((A_x/A6).^2 - 1/Mc_x.^2 * (2/(g+1) * (1 + (g-1)...
       * Mc_x.^2 /2)).^( (g+1)/(g-1)));
    Mc_x(i) = fzero(f,1);
end

%Supersonic Region

for i = find(Ac_x ==min(Ac_x)): length(Ac_x)
    A_x = Ac_x(i);
    f = @(Mc_x) ((A_x/A6).^2 - 1/Mc_x.^2 * (2/(g+1) * (1 + (g-1)...
       * Mc_x.^2 /2)).^( (g+1)/(g-1)));
    Mc_x(i) = fzero(f,10.5);
end

%Properties
Pc_x   = P05*(1+(g-1)*Mc_x.^2/2).^(-g/(g-1));
Tc_x   = T05*(1+(g-1)*Mc_x.^2/2).^(-1);
rhoc_x = rho05./((1+(g-1)/2*Mc_x.^2).^(1/(g-1)));
ac_x   = sqrt(g*R*Tc_x);
uc_x   = ac_x.*Mc_x;

%Plot

figure(7);
plot(xc,Rc_x,'k','LineWidth',2);
ylim([0 0.11]);
xlim ([0 0.258]);
title ('Conical Nozzle Geometry');
xlabel('Position (m)');
xlim([0 xc(end)]);
grid on

figure(8);
plot(xc,Mc_x,'k','LineWidth',2);
title('Mach Number - Conical Nozzle');
xlabel('Position (m)');
ylabel('Mach Number M');
xlim([0 xc(end)]);
grid on

figure(9);
plot(xc,Pc_x/Pc_x(1),'k','LineWidth',2);
title('Pressure Ratio - Conical Nozzle');
xlabel('Position (m)');
ylabel('Pressure ratio P/P0');
xlim([0 xc(end)]);
grid on

figure(10);
plot(xc,Tc_x/Tc_x(1),'k','LineWidth',2);
title('Temperature Ratio - Conical Nozzle');
ylabel('Temperature ratio T/T0');
xlim([0 xc(end)]);
grid on

figure(11);
plot(xc,rhoc_x,'k','LineWidth',2);
title('Density - Conical Nozzle');
ylabel('Density [kg/m^3]');
xlabel('Position (m)');
xlim([0 xc(end)]);
grid on;

figure(12);
plot(xc,uc_x,'k','LineWidth',2);
title('Flow Velocity - Conical Nozzle');
xlabel('Position (m)');
ylabel('Velocity (m/s)');
xlim([0 xc(end)]);
grid on