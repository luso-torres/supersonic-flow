function [xb,Rb_x,Mb_x,Pb_x,Tb_x,rhob_x,ub_x] = bell_geometry...
                                  (x_t,x_t2,Dt,P05,T05,rho05,A6)

g    = 1.41203;             %Degrees of freedom in a gas (gamma)            
R    = 287;                 %Gas Universal Constant [J*Kg^-1*K^-1]

%Auxiliar function to round numbers
roundn = @(x,n) round(x.*10.^n)./10.^n;


%Bell Shaped Geometry

x1 = 0:0.00001:0.09894;             %Geometric Restriction
x2 = x1(end):0.00001:0.152232;      %Geometric Restriction
x3 = x2(end):0.00001:roundn(x_t2,5);      
xb = [x1 x2 x3];

Rb_x1 =  sqrt(0.12^2-(x1-0).^2)-0.02;
Rb_x2 = -sqrt(Dt^2-(x2-x_t).^2)+2*Dt;
Rb_x3 =  sqrt(0.15^2-(x3-x_t2).^2)-0.05;
Rb_x  =  [Rb_x1 Rb_x2 Rb_x3];

Ab_x = Rb_x;
Mb_x = zeros(length(Ab_x),1);

%Mach at subsonic region

for i=1:find(Ab_x==min(Ab_x))
    A_x = Ab_x(i);
    f = @(Mb_x) ((A_x/A6).^2 - 1/Mb_x.^2 * (2/(g+1) * (1 + (g-1)...
       * Mb_x.^2 /2)).^( (g+1)/(g-1)));
    Mb_x(i) = fzero(f,1);
end

%Mach at supersonic region

for i =find(Ab_x == min(Ab_x)): length(Ab_x)
    A_x = Ab_x(i);
    f = @(Mb_x) ((A_x/A6).^2 - 1/Mb_x.^2 * (2/(g+1) * (1 + (g-1)...
       * Mb_x.^2 /2)).^( (g+1)/(g-1)));
    Mb_x(i) = fzero(f,10);
end

%Properties
Pb_x   = P05*(1+(g-1)*Mb_x.^2/2).^(-g/(g-1));
Tb_x   = T05*(1+(g-1)*Mb_x.^2/2).^(-1);
rhob_x = rho05./((1+((g-1)/2)*Mb_x.^2).^(1/(g-1)));
ab_x   = sqrt(g*R*Tb_x);
ub_x   = ab_x.*Mb_x;

figure(1);
plot(xb,Rb_x,'r','LineWidth',2);
ylim([0 0.11]);
title('Bell Shaped Nozzle Geometry');
xlabel('Position (m)');
xlim([0 xb(end)]);
grid on

figure(2);
plot(xb,Mb_x,'r','LineWidth',2);
title('Mach Number - Bell Shaped Nozzle');
xlabel('Position (m)');
ylabel('Mach number M');
xlim([0 xb(end)]);
grid on

figure(3);
plot(xb,Pb_x/Pb_x(1),'r','LineWidth',2);
title('Pressure Ratio - Bell Shaped Nozzle');
xlabel('Position (m)');
ylabel('Pressure ratio P/P0');
xlim([0 xb(end)]);
grid on

figure(4);
plot(xb, Tb_x/Tb_x(1),'r','LineWidth',2);
title('Temperature Ratio - Bell Shaped Nozzle');
ylabel('Temperature ratio T/T0');
xlabel('Position (m)');
xlim([0 xb(end)]);
grid on

figure(5);
plot(xb,rhob_x,'r','LineWidth',2);
title('Density - Bell Shaped Nozzle');
ylabel('Density [kg/m^3]');
xlabel('Position (m)');
xlim([0 xb(end)]);
grid on

figure(6);
plot(xb,ub_x,'r','LineWidth',2);
title('Flow Velocity - Bell Shaped Nozzle');
xlabel('Position (m)');
ylabel('Velocity (m/s)');
xlim([0 xb(end)]);
grid on

end