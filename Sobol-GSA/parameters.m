%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 15/09/2026                  %
% ----------------------------------------------------------- %
%  Main code: Sobol_main.m                                    %
%  Current function: parameters.m                             %
%    * parameter names & units, baseline values,              %
%      and defined lower and upper bounds                     %
%  Current model:                                             %
%      Bone fracture healing ODE model                        %
%      adapted from Trejo et al., 2019                        %
%      with chosen values (Lafuente-Gracia et al., 2026)      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p0,LB,UB,parNames,outputNames] = parameters

    % Parameter names
    parNames = {'ke1',      % Engulfing debris rate of M1 [/day]; Trejo range: [3 48]; value: 12
                'ke2',      % Engulfing debris rate of M2 [/day]; Trejo range: [3 48]; value: 24
                'aed',      % Half-saturation of debris [cells/mL]; value: 4.71e6
                'kmax',     % Maximal migration rate [/day]; Trejo range: [0.015 0.1]; value: 0.015
                'Mmax',     % Maximal macrophage density [cells/mL]; Trejo range: [6e5 1e6]; value: 1e6
                'k01',      % Activation rate of M0 to M1 [/day]; Trejo range: [0.55 0.611]; value: 0.611
                'k02',      % Activation rate of M0 to M2 [/day]; Trejo range: [0.0843 0.3]; value: 0.0836
                'a01',      % Half-saturation of c1 to activate M1 [ng/mL]; value: 0.01
                'a02',      % Half-saturation of c2 to activate M2 [ng/mL]; value: 0.005
                'k12',      % Transition rate from M1 to M2 [/day]; Trejo range: [0.075 0.083]; value: 0.075
                'k21',      % Transition rate from M2 to M1 [/day]; Trejo range: [0.005 0.05]; value: 0.05
                'd0',       % Apoptosis rate of M0 [/day]; Trejo range: [0.156 0.2]; value: 0.156
                'd1',       % Apoptosis rate of M1 [/day]; Trejo range: [0.121 0.2]; value: 0.121
                'd2',       % Apoptosis rate of M2 [/day]; Trejo range: [0.163 0.2]; value: 0.163
                'k0',       % Secretion rate of c1 by debris [ng/(cells*day)]; Trejo range: [5e-7 8.5e-6]; value: 5e-7
                'k1',       % Secretion rate of c1 by M1 [ng/(cells*day)]; value: 8.3e-6
                'k2',       % Secretion rate of c2 by M2 [ng/(cells*day)]; value: 3.72e-6
                'k3',       % Secretion rate of c2 by Cs [ng/(cells*day)]; Trejo range: [7e-7 8e-6]; value: 8e-6
                'dc1',      % Decay rate of c1 [/day]; Trejo range: [12.79 55]; value: 12.79
                'dc2',      % Decay rate of c2 [/day]; Trejo range: [2.5 4.632]; value: 4.632
                'a12',      % Effectiveness of c2 inhibition of c1 synthesis [ng/mL]; value: 0.025
                'a22',      % Effectiveness of c2 inhibition of c2 synthesis [ng/mL]; value: 0.1
                'aps',      % Effectiveness of c1 inhibition of Cs proliferation [ng/mL]; value: 3.162
                'asb1',     % Effectiveness of c1 inhibition of Cs differentiation [ng/mL]; value: 0.1
                'apb',      % Effectiveness of c1 inhibition of Cb proliferation [ng/mL]; value: 10
                'aps1',     % Constant enhancement of c1 to Cs proliferation [ng/mL]; value: 20
                'kps',      % Cs proliferation rate [/day]; value: 0.5
                'kls',      % Cs carrying capacity [cells/mL]; value: 1e6
                'ds',       % Differentiation rate of Cs into Cb [/day]; value: 1
                'kpb',      % Cb proliferation rate [/day]; value: 0.2202
                'klb',      % Cb carrying capacity [cells/mL]; value: 1e6
                'db',       % Differentiation rate of Cb into osteocytes [/day]; value: 0.15
                'pcs',      % Fibrocartilage synthesis rate [g/(cells*day)]; value: 3e-6
                'qcd1',     % Fibrocartilage degradation rate [mL/(cells*day)]; value: 3e-6
                'qcd2',     % Fibrocartilage degradation rate by osteoclasts [mL/(cells*day)]; value: 0.2e-6
                'pbs',      % Woven bone synthesis rate [g/(cells*day)]; value: 5e-8
                'qbd'};     % Woven bone degradation rate [mL/(cells*day)]; value: 5e-8
    
    % Baseline values
    p0 = [12;           % ke1
         24;            % ke2
         4.71e6;        % aed
         0.015;         % kmax
         1e6;           % Mmax
         0.611;         % k01
         0.0836;        % k02
         0.01;          % a01
         0.005;         % a02
         0.075;         % k12
         0.05;          % k21
         0.156;         % d0
         0.121;         % d1
         0.163;         % d2
         5e-7;          % k0
         8.3e-6;        % k1
         3.72e-6;       % k2
         8e-6;          % k3
         12.79;         % dc1
         4.632;         % dc2
         0.025;         % a12
         0.1;           % a22
         3.162;         % aps
         0.1;           % asb1
         10;            % apb
         20;            % aps1
         0.5;           % kps
         1e6;           % kls
         1;             % ds
         0.2202;        % kpb
         1e6;           % klb
         0.15;          % db
         3e-6;          % pcs
         3e-6;          % qcd1
         0.2e-6;        % qcd2
         5e-8;          % pbs
         5e-8];         % qbd
    
    % Lower and upper bounds with +-50%
    LB = 0.5*p0;
    UB = 1.5*p0;
    
    % Sanity check
    ParamTable = table(parNames,p0,LB,UB,'VariableNames',{'Parameter','Nominal','LowerBound','UpperBound'});
    writetable(ParamTable,'output/ParameterRanges.xlsx');
    
    % Outputs of the sensitivity analysis
    outputNames = {'M1 peak';
                   'M2 peak';
                   'c1 peak';
                   'c2 peak';
                   'Cs peak';
                   'Cb peak';
                   'mf peak';
                   'mb peak';
                   'M1 time-to-peak';
                   'M2 time-to-peak';
                   'c1 time-to-peak';
                   'c2 time-to-peak';
                   'Cs time-to-peak';
                   'Cb time-to-peak';
                   'mf time-to-peak';
                   %'mb time-to-peak'
                   };

end