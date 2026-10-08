%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 15/09/2026                  %
% ----------------------------------------------------------- %
%  Main code: Sobol_main.m                                    %
%  Current function: evaluate_model.m                         %
%    * defines initial conditions                             %
%              & ODEs (see function at the bottom)            %
%    * solves ODE model                                       %
%    * saves outputs                                          %
%    * saves outputs                                          %
%  Current model:                                             %
%      Bone fracture healing ODE model                        %
%      adapted from Trejo et al., 2019                        %
%      with chosen values (Lafuente-Gracia et al., 2026)      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out,outVec] = evaluate_model(p)

    time = [0 400];
    
    % Initial conditions
    y0 = zeros(10,1);
    y0(1) = 5e7;              % Initial D density (= density of debis/necrotic cells) [cells/mL]
                              %%% Range of values (Trejo): DIC = [1*10^6 2*10^8];
                              %%% section 6.4 - p.12: D(0)=3*10^5 simple, D(0)=5*10^7 moderate, & D(0)=2*10^8 severe fracture
                              %%% section 6.5.1 - p.13: D(0)=3*10^5, D(0)=2*10^7, D(0)=5*10^7
    y0(2) = 4000;             % Initial M0 density (unactivated macrophages) [cells/mL]
    y0(3) = 0;                % Initial M1 density (classical macrophages) [cells/mL]
    y0(4) = 0;                % Initial M2 density (alternative macrophages) [cells/mL]
    y0(5) = 1;                % Initial c1 concentration (pro-inflammatory cytokines) [ng/mL]
    y0(6) = 0;                % Initial c2 concentration (anti-inflammatory cytokines) [ng/mL]
                              %%% section 6.5.1 - p.13: c2(0)=0, c2(0)=10, c2(0)=100
    y0(7) = 1000;             % Initial Cs density (SSPCs) [cells/mL]
    y0(8) = 0;                % Initial Cb density (osteoblasts) [cells/mL]
    y0(9) = 0;                % Initial mc density (fibrocartilage) [g/mL]
    y0(10) = 0;               % Initial mb density (bone) [g/mL]
    
    % Solve system
    [T,Y] = ode23s(@(t,y) odefun(t,y,p),time,y0);
    
    % Peak values
    [M1_peak,M1_idx] = max(Y(:,3));
    [M2_peak,M2_idx] = max(Y(:,4));
    [c1_peak,c1_idx] = max(Y(:,5));
    [c2_peak,c2_idx] = max(Y(:,6));
    [Cs_peak,Cs_idx] = max(Y(:,7));
    [Cb_peak,Cb_idx] = max(Y(:,8));
    [mf_peak,mf_idx] = max(Y(:,9));
    [mb_peak,mb_idx] = max(Y(:,10));
    
    % Time to peak
    M1_tpeak = T(M1_idx);
    M2_tpeak = T(M2_idx);
    c1_tpeak = T(c1_idx);
    c2_tpeak = T(c2_idx);
    Cs_tpeak = T(Cs_idx);
    Cb_tpeak = T(Cb_idx);
    mf_tpeak = T(mf_idx);
    %mb_tpeak = T(mb_idx);

    % Output structure
    out.M1_peak = M1_peak;
    out.M2_peak = M2_peak;
    out.c1_peak = c1_peak;
    out.c2_peak = c2_peak;
    out.Cs_peak = Cs_peak;
    out.Cb_peak = Cb_peak;
    out.mf_peak = mf_peak;
    out.mb_peak = mb_peak;
    out.M1_tpeak = M1_tpeak;
    out.M2_tpeak = M2_tpeak;
    out.c1_tpeak = c1_tpeak;
    out.c2_tpeak = c2_tpeak;
    out.Cs_tpeak = Cs_tpeak;
    out.Cb_tpeak = Cb_tpeak;
    out.mf_tpeak = mf_tpeak;
    %out.mb_tpeak = mb_tpeak;

    outVec = [ ...
        out.M1_peak ...
        out.M2_peak ...
        out.c1_peak ...
        out.c2_peak ...
        out.Cs_peak ...
        out.Cb_peak ...
        out.mf_peak ...
        out.mb_peak ...
        out.M1_tpeak ...
        out.M2_tpeak ...
        out.c1_tpeak ...
        out.c2_tpeak ...
        out.Cs_tpeak ...
        out.Cb_tpeak ...
        out.mf_tpeak ...
        %out.mb_tpeak
        ];

end

function dydt = odefun(~,y,p)

    % unpack parameters
    ke1  = p(1);
    ke2  = p(2);
    aed  = p(3);
    kmax = p(4);
    Mmax = p(5);
    k01 = p(6);
    k02 = p(7);
    a01 = p(8);
    a02 = p(9);
    k12 = p(10);
    k21 = p(11);
    d0 = p(12);
    d1 = p(13);
    d2 = p(14);
    k0 = p(15);
    k1 = p(16);
    k2 = p(17);
    k3 = p(18);
    dc1 = p(19);
    dc2 = p(20);
    a12 = p(21);
    a22 = p(22);
    aps  = p(23);
    asb1 = p(24);
    apb  = p(25);
    aps1 = p(26);
    kps = p(27);
    kls = p(28);
    ds = p(29);
    kpb = p(30);
    klb = p(31);
    db = p(32);
    pcs  = p(33);
    qcd1 = p(34);
    qcd2 = p(35);
    pbs = p(36);
    qbd = p(37);

    % terms
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

    % ODEs
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