%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Bone fracture healing: inflammation ODE model   %
%           adapted from Trejo et al., 2019         %
%        Implemented by Laura Lafuente-Gracia       %
%             Last revision: 23/07/2025             %
%           Male vs. female parameter sets          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

%% Parameter values with descriptions and units
E = 'E1';                           % Choose equilibria: E0 (non-healing), E1 (healing), E2 (non-/delayed-healing)
                                    % E1 = reference state
maxinf = 50;
maxrep = 400;

time = [0 maxrep];                  % Total simulated time [days]
                                    %%% Range of values (Trejo): time = [0 400]

% New parameters male
ke1 = 12;                           % Engulfing debris rate of M1 [/day]
                                    %%% Range of values (Trejo): [3 48]
ke2 = 24;                           % Engulfing debris rate of M2 [/day]
                                    %%% Range of values (Trejo): [3 48]
aed = 4.71*10^6;                    % Half-saturation of debris [cells/mL]
kmax = 0.015;                       % Maximal migration rate [/day]
                                    %%% Range of values (Trejo): [0.015 0.1]  
Mmax = 1*10^6;                      % Maximal macrophages density [cells/mL]
                                    %%% Range of values (Trejo): [6*10^5 1*10^6]
k01 = 0.611;                        % Activation rate of M1 [/day]
                                    %%% Range of values (Trejo): [0.55 0.611]
k02 = 0.0836;                       % Activation rate of M0 to M2 [/day]
                                    %%% Range of values (Trejo): [0.0843 0.3]
a01 = 0.01;                         % Half-saturation of c1 to activate M1 [ng/mL] 
a02 = 0.005;                        % Half-saturation of c2 to activate M2 [ng/mL]
k12 = 0.075;            	        % Transition rate from M1 to M2 [/day]
                                    %%% Range of values (Trejo): [0.075 0.083]
k21 = 0.05;                         % Transition rate from M2 to M1 [/day]
                                    %%% Range of values (Trejo): [0.005 0.05]
d0 = 0.156;                         % Emigration rate of M0 [/day]
                                    %%% Range of values (Trejo): [0.156 0.2]
d1 = 0.121;                         % Emigration rate of M1 [/day]
                                    %%% Range of values (Trejo): [0.121 0.2]
d2 = 0.163;                         % Emigration rate of M2 [/day]
                                    %%% Range of values (Trejo): [0.163 0.2]
k0 = 5*10^(-7);                     % Secretion rate of c1 by debris [ng/cells/day]
                                    %%% Range of values (Trejo): [5*10^(-7) 8.5*10^(-6)]
k1 = 8.3*10^(-6);                   % Secretion rate of c1 by M1 macrophages  [ng/cells/day] 
k2 = 3.72*10^(-6);                  % Secretion rate of c2 by M2 macrophages  [ng/cells/day]
k3 = 8*10^(-6);                     % Secretion rate of c2 by SSPCs [ng/cells/day]
                                    %%% Range of values (Trejo): [7*10^(-7) 8*10^(-6)]
dc1 = 12.79;                        % Decay rate of c1 [/day]
                                    %%% Range of values (Trejo): [12.79 55]
dc2 = 4.632;                        % Decay rate of c2 [/day]
                                    %%% Range of values (Trejo): [2.5 4.632]
a12 = 0.025;                        % Effectiveness of c2 inhibition of c1 synthesis[ng/mL]
a22 = 0.1;                          % Effectiveness of c2 inhibition of c2 synthesis [ng/mL] 
aps = 3.162;                        % Effectiveness of c1 inhibition of Cs proliferation [ng/mL]
asb1 = 0.1;                         % Effectiveness of c1 inhibition of Cs differentiation [ng/mL] 
apb = 10;                           % Effectiveness of c1 inhibition of Cb proliferation  ng/mL 
aps1 = 20;                          % Constant enhancement of c1 to Cs proliferation  ng/mL 
kps = 0.5;                          % Proliferation rate of Cs [/day]
kls = 1*10^6;                       % Carrying capacity of Cs [cells/mL]
ds = 1;                             % Differentiation rate of Cs [/day] 
%ds =1 for E0, =1 for E1, =0.1 for E2
kpb = 0.2202;                       % Proliferation rate of Cb [/day]
klb = 1*10^6;                       % Carrying capacity of Cb [cells/mL]
db = 0.15;                          % Differentiation rate of Cb [/day] to osteocytes + apoptosis
%db =0.3 for E0, =0.15 for E1, =0.15 for E2
pcs = 3*10^(-6);                    % Fibrocartilage synthesis rate [g/cells/day]
qcd1 = 3*10^(-6);                   % Fibrocartilage degradation rate [mL/cells/day]
qcd2 = 0.2*10^(-6);                 % Fibrocartilage degradation rate by osteoclasts [mL/cells/day]
pbs = 5*10^(-8);                    % Bone tissue synthesis rate [g/cells/day]
qbd = 5*10^(-8);                    % Bone tissue degradation rate [mL/cells/day]

switch E
    case 'E0'
        disp('*** Equilibrium point E0: non-healing ***');
        ds = 1;     % E0 = E1
        db = 0.3;   % E0 = E1*2
    case 'E1'
        disp('*** Equilibrium point E1: healing ***');
        ds = 1;
        db = 0.15;
    case 'E2'
        disp('*** Equilibrium point E2: non- or delayed-healing ***');
        ds = 0.1;   % E2 = E1/10
        db = 0.15;  % E2 = E1
end


%% Initial conditions

y0 = zeros(10,1);
y0(1) = 5*10^7;                     % Initial D density (= density of debis/necrotic cells) [cells/mL]
                                    %%% Range of values (Trejo): DIC = [1*10^6 2*10^8];
                                    %%% section 6.4 - p.12: D(0)=3*10^5 simple, D(0)=5*10^7 moderate, & D(0)=2*10^8 severe fracture
                                    %%% section 6.5.1 - p.13: D(0)=3*10^5, D(0)=2*10^7, D(0)=5*10^7
y0(2) = 4000;                       % Initial M0 density (unactivated macrophages) [cells/mL]
y0(3) = 0;                          % Initial M1 density (classical macrophages) [cells/mL]
y0(4) = 0;                          % Initial M2 density (alternative macrophages) [cells/mL]
y0(5) = 1;                          % Initial c1 concentration (pro-inflammatory cytokines) [ng/mL]
y0(6) = 0;                          % Initial c2 concentration (anti-inflammatory cytokines) [ng/mL]
                                    %%% section 6.5.1 - p.13: c2(0)=0, c2(0)=10, c2(0)=100
y0(7) = 1000;                       % Initial Cs density (SSPCs) [cells/mL]
y0(8) = 0;                          % Initial Cb density (osteoblasts) [cells/mL]
y0(9) = 0;                          % Initial mc density (fibrocartilage) [g/mL]
y0(10) = 0;                         % Initial mb density (bone) [g/mL]


%% Solve ODE system

% Base model
[T,Y] = ode23s(@(t,y) odefun(t,y,ke1,ke2,aed,kmax,Mmax,k01,k02,k12,k21,d0,d1,d2,k0,k1,k2,k3,dc1,dc2,a12,a22,aps,asb1,a01,a02,apb,aps1,kps,ds,kpb,db,pcs,qcd1,qcd2,pbs,qbd,klb,kls), time, y0);

As = kps * (aps^2 + aps1*Y(:,5)) ./ (aps^2 + Y(:,5).^2);
F1 = ds * asb1 ./ (asb1 + Y(:,5));
Ab = kpb * apb/(apb + Y(:,5));

% New parameters female
% inflammatory: +-50%
Mmax = Mmax*1.5;
k01 = k01*1.5;
k02 = k02*1.5;
ke1 = ke1*1.5;
ke2 = ke2*1.5;
k1 = k1*0.5;
k2 = k2*1.5;
k3 = k3*1.5;
% repair: +-25%
kls = kls*0.75;
kps = kps*1.25;
ds = ds*0.75;
klb = klb*0.75;
kpb = kpb*1.25;

[Tnew,Ynew] = ode23s(@(tnew,ynew) odefun(tnew,ynew,ke1,ke2,aed,kmax,Mmax,k01,k02,k12,k21,d0,d1,d2,k0,k1,k2,k3,dc1,dc2,a12,a22,aps,asb1,a01,a02,apb,aps1,kps,ds,kpb,db,pcs,qcd1,qcd2,pbs,qbd,klb,kls), time, y0);

As_new = kps * (aps^2 + aps1*Ynew(:,5)) ./ (aps^2 + Ynew(:,5).^2);
F1_new = ds * asb1 ./ (asb1 + Ynew(:,5));
Ab_new = kpb * apb/(apb + Ynew(:,5));


%% Plot solution & save images

%%% Inflammatory response
Inflam = figure;
ir = uipanel('Parent',Inflam,'BorderType','none'); 
ir.Title = 'Inflammatory response';
ir.TitlePosition = 'centertop'; 
ir.FontSize = 12;
ir.FontWeight = 'bold';
Inflam.Position = [10 100 900 400];

subplot(2,3,1,'Parent',ir,'Color', 'white')
plot(T,Y(:,1), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,1), 'LineWidth', 1);
title('Debris');
xlabel('time [days]');
ylabel('density [cells/mL]');
legend('male', 'female');
xlim([0 maxinf]);

subplot(2,3,2,'Parent',ir)
plot(T,Y(:,5), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,5), 'LineWidth', 1);
title('c1 cytokines');
xlabel('time [days]');
ylabel('concentration [ng/mL]');
xlim([0 maxinf]);

subplot(2,3,3,'Parent',ir)
plot(T,Y(:,6), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,6), 'LineWidth', 1); 
title('c2 cytokines');
xlabel('time [days]');
ylabel('concentration [ng/mL]');
xlim([0 maxinf]);

subplot(2,3,4,'Parent',ir)
plot(T,Y(:,2), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,2), 'LineWidth', 1); 
title('M0 macrophages');
xlabel('time [days]');
ylabel('density [cells/mL]');
xlim([0 maxinf]);

subplot(2,3,5,'Parent',ir)
plot(T,Y(:,3), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,3), 'LineWidth', 1);
title('M1 macrophages');
xlabel('time [days]');
ylabel('density [cells/mL]');
xlim([0 maxinf]);

subplot(2,3,6,'Parent',ir)
plot(T,Y(:,4), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,4), 'LineWidth', 1); 
title('M2 macrophages');
xlabel('time [days]');
ylabel('density [cells/mL]');
xlim([0 maxinf]);

saveas(gcf,'inflammatory-response','png')

%%% Repair phase
Repair = figure;
rp = uipanel('Parent',Repair,'BorderType','none');
rp.Title = 'Repair phase';
rp.TitlePosition = 'centertop'; 
rp.FontSize = 12;
rp.FontWeight = 'bold';

subplot(2,2,1,'Parent',rp)
plot(T,Y(:,7), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,7), 'LineWidth', 1); 
title('SSPCs');
xlabel('time [days]');
ylabel('density [cells/mL]');
legend('male', 'female');
xlim([0 100]);

subplot(2,2,2,'Parent',rp)
plot(T,Y(:,8), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,8), 'LineWidth', 1);
title('Osteoblasts');
xlabel('time [days]');
ylabel('density [cells/mL]');
xlim([0 100]);

subplot(2,2,3,'Parent',rp)
plot(T,Y(:,9), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,9), 'LineWidth', 1); 
title('Fibrocartilage');
xlabel('time [days]');
ylabel('density [g/mL]');
xlim([0 100]);

subplot(2,2,4,'Parent',rp)
plot(T,Y(:,10), 'LineWidth', 1);
hold on
plot(Tnew,Ynew(:,10), 'LineWidth', 1);
title('Bone');
xlabel('time [days]');
ylabel('density [g/mL]');
xlim([0 maxrep]);

saveas(gcf,'repair-phase','png')


%% ODE system: definition of equations & their terms

function dydt = odefun(~,y,ke1,ke2,aed,kmax,Mmax,k01,k02,k12,k21,d0,d1,d2,k0,k1,k2,k3,dc1,dc2,a12,a22,aps,asb1,a01,a02,apb,aps1,kps,ds,kpb,db,pcs,qcd1,qcd2,pbs,qbd,klb,kls)

    RD = y(1)/(aed+y(1));                       % Debris engulfing rate
    M  = y(2) + y(3) + y(4);                    % Total density of macrophages
    RM = kmax * (1 - M/Mmax) * y(1);            % Migration rate of unactivated macrophages
    G1 = k01 * y(5)/(a01 + y(5));               % Differentiation rate of M1
    G2 = k02 * y(6)/(a02 + y(6));               % Differentiation rate of M2
    H1 = a12/(a12+y(6));                        % Inhibition of c1
    H2 = a22/(a22+y(6));                        % Inhibition of c2
    As = kps * (aps^2 + aps1*y(5)) / (aps^2 + y(5)^2);  % Proliferation of SSPCs
    F1 = ds * asb1/(asb1 + y(5));               % Differentiation of SSPCs to osteoblasts
    Ab = kpb * apb/(apb + y(5));                % Proliferation of osteoblasts

    dydt = zeros(10,1);
    dydt(1)  = -RD * ( ke1*y(3) + ke2*y(4) );                   % Debris (D)
    dydt(2)  = RM - G1*y(2) - G2*y(2) - d0*y(2);                % Unactivated macrophages (M0)
    dydt(3)  = G1*y(2) + k21*y(4) - k12*y(3) - d1*y(3);         % Classical macrophages (M1)
    dydt(4)  = G2*y(2) + k12*y(3) - k21*y(4) - d2*y(4);         % Alternative macrophages (M1)      
    dydt(5)  = H1 * ( k0*y(1) + k1*y(3) ) - dc1*y(5);           % Pro-inflammatory cytokines (c1)
    dydt(6)  = H2 * ( k2*y(4) + k3*y(7) ) - dc2*y(6);           % Anti-inflammatory cytokines (c2)
    dydt(7)  = As*y(7) * ( 1 - y(7)/kls ) - F1*y(7);            % SSPCs (Cs)
    dydt(8)  = Ab*y(8) * ( 1 - y(8)/klb ) + F1*y(7) - db*y(8);  % Osteoblasts (Cb)
    dydt(9)  = (pcs - qcd1*y(9)) * y(7) - qcd2*y(9)*y(8);       % Fibrocartilage (mc)
    dydt(10) = (pbs - qbd*y(10)) * y(8);                        % Bone (mb)
end