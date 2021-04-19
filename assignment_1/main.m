% If you want to execute one exercise separately, press ctrl+enter in the
% desired cell

% Remember to save the script before executing it!

%% Exercise 1 - Simulation stability

clc                             % Clean command window
clear                           % Everything in the workspace is deleted

% Create a time vector
dt=2.51e-3;                         % The time step is set to 0.1 s
T=1;                            % The stop time is set to 1 s
t=[0:dt:T];                     % A time vector is created with start time = 0, 
                                % time step = dt and stop time = T

% Create a for loop
x(1)=-5;                         % Sets a start value for the x-vector

% step 1
Kc = 8e-10;
Be = 1000e6;
V = 1e-3;
p_lam(1) = 100e3;
p0 = 10e6;

%Laminar flow
for i=1:(T/dt)
    q_in_lam(i) = Kc*(p0-p_lam(i));
    p_lam(i+1) = p_lam(i)+q_in_lam(i)*Be/V*dt;
end                             

p_turb(1) = 100e3;
Cq = 0.67;
d = 2e-3;
A = d^2/4*pi;
rho = 860;
%Turbulent flow
for i=1:(T/dt)
    if p0-p_turb(i) > 0
        q_in_turb(i) = Cq*A*sqrt((p0-p_turb(i))*2/rho);
    else
        q_in_turb(i) = -Cq*A*sqrt((p_turb(i)-p0)*2/rho);
    end
    p_turb(i+1) = p_turb(i)+q_in_turb(i)*Be/V*dt;
end

figure()
plot(t,p_lam)                       % The x-vector is plotted a against the time vector
xlabel('Time [s]')                  % Defines a name to the x-axis
ylabel('Pressure [Pa]') % Defines a name to the y-axis
title('Laminar flow')

figure()
plot(t,p_turb)                       % The x-vector is plotted a against the time vector
xlabel('Time [s]')                  % Defines a name to the x-axis
ylabel('Pressure [Pa]')
title('Turbulent flow')


%% Exercise 2 - Time domain analysis

% Execute this cell and open the Simulink file. You can change values in
% the Simulink file from here.

%% Plotting exercise 2

%Data
clear
Ps = 20e6;
xv = 1e-3;
w = 10e-3;
Cq = 0.67;
A1 = 50e-4;
A2 = A1;
Vt = 50e-3;
Be = 1000e6;
rho = 860;
M = 1000;
F = 50e3;
Bp = 10e3;

% This cell cannot be executed before simulating "Simulink_example"
modelPath = 'Simulink_example';

% Either open the simulink file and simulate:
open(modelPath)

% Or simulate simulink model from script:
Stop_Time = 5;
sim(modelPath,Stop_Time)

% Plotting the flow from the block "simout" in Simulink_example
figure()
plot(Results.time,Results.signals.values(:,1),'k','LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Pressure [Pa]')
hold on                 % The old plot will be kept if you change parmeter values
                        % and want to plot again

% Plotting the pressure from the block "simout" in Simulink_example
figure()
plot(Results.time,Results.signals.values(:,2),'k','LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Flow [m^3/s]')
%ylim([0 1.1e7])             %Specifying the upper and lower limit of the y-axis
hold on

figure()
plot(Results.time,Results.signals.values(:,3),'k','LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Xp_{dot} [m]')
%ylim([0 1.1e7])             %Specifying the upper and lower limit of the y-axis
hold on


%% Exercise 3 - Frequency domain analysis 
%Variables
s = tf('s');
Ps = 15e6;
d = 1e-3;
Cq = 0.67;
V = 1e-3;
Be = 1000e6;
rho = 860;
Dm = 20e-6;
Jm = 0.02;
Bm = 0.13;
Te = 20;
PL = 7.2e6;

kq = Cq*2*d*pi/4*sqrt(2/rho*(Ps-PL));
kc = Cq*d^2*pi/4/(2*sqrt(rho/2*(Ps-PL)));

konst = 2*pi*kc*Bm+Dm^2/(2*pi);
w_h = sqrt(Be/(2*pi*V*Jm)*konst);
delta_h = (pi*V*Bm/Be+pi*kc*Jm)*w_h/konst;
K = Dm/konst;

d = 0.1e-3;
Te = 1;
Nm = (kq*d-Te*(2*pi*(V*s/Be+kc)/Dm))*K/((s/w_h)^2+2*delta_h/w_h*s+1);

Nm_eval = freqresp(Nm,0);

d = 1e-3;
Te = 20;

%With Te=0
G0 = kq*K/((s/w_h)^2+2*delta_h/w_h*s+1);
figure()
bode(G0)

%With d=0
Gd = (2*pi*(V*s/Be+kc)/Dm)*K/((s/w_h)^2+2*delta_h/w_h*s+1);
figure()
bode(Gd)

wb = bandwidth(G0);
%wb = 1000;
sim('mainSim')

figure()
plot(Nm_sim)
legend('Non-lin Simulink','Linear Analytical')
ylabel('N_m [rad/s]')
title('Motor speed')
