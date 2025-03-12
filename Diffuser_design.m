% Atmospheric conditions @ 22km
clc; close all
[T, a, P, rho] = atmosisa(22000);
M1 = 3.5;     %Mach Number

%Free-stream conditions

R = 287;      %Air gas constant [J kg^-1 K^-1]
P1 = P;       %Free-stream pressure [Pa = N/m^2]
T1 = T;       %Free-stream temperature [K]
rho1 = rho;   %Free-stream density [kg/m^3]
a1 = a;       %Speed of sound [m/s]
V1 = M1*a1;   %Free-stream velocity [m/s]
g = 1.41203;  %Adiabatic coefficient

P01 = P1*((1+(g-1)/2*M1^2)^(g/(g-1))); %Total pressure
cp = (g*R)/(g-1); %1004.5 [J kg^-1 K^-1];

Pm = 0;        %loop variable

l1 = .111994;
l2 = .029684;

for beta1 = 1:.1:89 %accuracy criteria
    Mn1  = M1*sind(beta1);% sin (degrees)
    Mn2  = sqrt((1+((g-1)/2)*(Mn1^2))/(g*Mn1^2-((g-1)/2)));
    rho2 = rho1*((g+1)*(Mn1^2))/(2+(g-1)*Mn1^2);
    theta1 = atand((2*cotd(beta1)*(((M1^2)*(sind(beta1))^2)-1))/...
             ((M1^2)*(g+cosd(2*beta1))+2)); %atand -> degrees;
    M2  = Mn2/sind(beta1-theta1);
    P2  = P1*(1+(2*g/(g+1))*(Mn1^2-1));
    T2  = T1*(P2/P1*rho1/rho2);
    ds  = cp*log(T2/T1)-R*log(P2/P1); %log -> ln function;
    P02 = P01*exp(-ds/R);
    a2  = sqrt(g*R*T2);
    V2  = M2*a2;


    for beta2 = 1:.1:89 %accuracy criteria;
        Mn2  = M2*sind(beta2);% sin (degrees)
        Mn3  = sqrt((1+((g-1)/2)*(Mn2^2))/(g*(Mn2^2)-...
        ((g-1)/2)));
        rho3  = rho2*((g+1)*(Mn2^2))/(2+(g-1)*Mn2^2);
        theta2 = atand((2*cotd(beta2)*(((M2^2)*sind(beta2)^2)-1))/...
                 ((M2^2)*(g+cosd(2*beta2))+2)); %atand -> degrees;
        M3  = Mn3/sind(beta2-theta2);
        P3  = P2*(1+(2*g/(g+1))*(Mn2^2-1));
        T3  = T2*(P3/P2*rho2/rho3);
        ds  = cp*log(T3/T2)-R*log(P3/P2); %log -> ln function;
        P03 = P02*exp(-ds/R);
        a3  = sqrt(g*R*T3);
        V3  = M3*a3;
    

    beta3 = 90;

    Mn3  = M3*sind(beta3);% sin (degrees)
    Mn4  = sqrt((1+((g-1)/2)*(Mn3^2))/(g*(Mn3^2)-((g-1)/2)));
    rho4 = rho3*((g+1)*(Mn3^2))/(2+(g-1)*Mn3^2);
    M4   = Mn4;
    P4   = P3*(1+(2*g/(g+1))*(Mn3^2-1));
    T4   = T3*(P4/P3*rho3/rho4);
    ds   = cp*log(T4/T3)-R*log(P4/P3); %log -> ln function;
    P04  = P03*exp(-ds/R);
    a4   = sqrt(g*R*T4);
    V4   = M4*a4;    
    
% Constraint
    if (Pm < P04 && Mn1>=1 && Mn2>=1 && M3>=1)
        Pm = P04;
        opt_M = [M1 M2 M3 M4];
        opt_P = [P1 P2 P3 P4];
        opt_T = [T1 T2 T3 T4];
        opt_a = [a1 a2 a3 a4];
        opt_V = [V1 V2 V3 V4];
        opt_rho = [rho1 rho2 rho3 rho4];
        opt_P0 = [P01 P02 P03 P04];
        opt_theta = [theta1 theta2];
        opt_beta = [beta1 beta2 beta3];
        n = P04/P01;
        D = 2*l1*sind(theta1)*(P2-P1)+2*l2*sind(theta2)*(P3-P2);
    end
    end    
end

figure(1)
x1 = 0:3;
plot(x1,opt_M)
hold on;
plot(x1,opt_M,'rd');

title ('Mach number Variation');
xlabel('Station');
ylabel('Mach Number');
set(gca, 'XTick', 0:4)

figure(2)
plot(x1,opt_P)
hold on
plot(x1,opt_P,'rd');

title ('Pressure Variation');
xlabel('Station');
ylabel('Pressure [Pa]');
set(gca, 'XTick', 0:3)


figure(3)
plot(x1,opt_T)
hold on
plot(x1,opt_T,'rd');

title ('Temperature Variaton');
xlabel('Station');
ylabel('Temperature [K]');
set(gca, 'XTick', 0:3)
        
figure(4)
plot(x1,opt_V)
hold on
plot(x1,opt_V,'rd');

title ('Velocity Variation');
xlabel('Station');
ylabel('Velocity [m/s]');
set(gca, 'XTick', 0:3)

figure(5)
plot(x1,opt_rho)
hold on
plot(x1,opt_rho,'rd');

title ('Density Variation');
xlabel('Station');
ylabel('Density [kg/m^3]');
set(gca, 'XTick', 0:3)