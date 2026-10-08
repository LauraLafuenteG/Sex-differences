%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 21/09/2026                  %
% ----------------------------------------------------------- %
%  Functions needed:                                          %
%    * evaluate_model.m: ODE model definition + ICs           %
%    * parameters.m: ODE parameter values & ranges for SA     %
%    * sobol_sampling.m                                       %
%    * sobol_indices.m                                        %
%    * plot_sobol.m                                           %
% Additional main code:                                       %
%    * Sobol_convergence.m: check convergence to choose N     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

if ~exist('output','dir')
    mkdir('output')
end

tic 

% Load model parameters
[p0,LB,UB,parNames,outputNames] = parameters;

% Sobol settings
D = length(p0);   % number of parameters
nOutputs = length(outputNames);
N = 8192;   % number of samples
            % power of 2 to preserve its low-discrepancy balance properties
nSim = N*(D+2);
fprintf('Number of simulations: %d\n',nSim);

% Generate samples
[A,B] = sobol_sampling(N,LB,UB);   % dimension: N x D (i.e. each row = one complete parameter set)

% Storage
YA  = zeros(N,nOutputs);      % model output obtained from A
YB  = zeros(N,nOutputs);      % model output obtained from B
YAB = zeros(N,nOutputs,D);    % model output obtained from mixed matrix AB

% Evaluate A
for i = 1:N
    [~,outVec] = evaluate_model(A(i,:));
    YA(i,:) = outVec;
end

% Evaluate B
for i = 1:N
    [~,outVec] = evaluate_model(B(i,:));
    YB(i,:) = outVec;
end

% Evaluate mixed matrices
% For every parameter k, we construct the mixed matrix AB^k: 
% all columns from A, except parameter k which is copied from B
for k = 1:D

    fprintf('Parameter %d/%d\n',k,D)

    AB = A;
    AB(:,k) = B(:,k);

    for i = 1:N
        [~,outVec] = evaluate_model(AB(i,:));
        YAB(i,:,k) = outVec;
    end

end

% Compute Sobol indices
Si  = zeros(D,nOutputs);   % main-effect index (first-order index, direct contribution)
                           % i.e. Si = fraction of output variance explained by parameter i alone
                           % e.g. Si = 0.5 --> 50% of the output variance is explained directly by that parameter

STi = zeros(D,nOutputs);   % total-effect index (total-order index, total influence)
                           % it includes direct effect, interactions, and nonlinear effects involving that parameter

% Example: Si = 0.2, STi = 0.65
% 65% total influence, of which 20% variance: parameter alone
%                               65% - 20% = 45% variance through higher-order interactions
% thus (0.65-0.2)/0.65 = 0.69 --> 69% of the parameter influence comes from interactions

for outID = 1:nOutputs
    [Si(:,outID),STi(:,outID)] = sobol_indices(YA(:,outID), YB(:,outID), squeeze(YAB(:,outID,:)));
end

save('output/Sobol_results.mat','Si','STi','parNames','outputNames','N');


%% Export data & plots per output
filename = sprintf('output/Sobol_summary_N%d.xlsx',N);

for outID = 1:nOutputs

    Interaction = STi(:,outID) - Si(:,outID);
    InteractionFraction = Interaction ./ max(STi(:,outID),1e-12);

    T = table(parNames(:), Si(:,outID), STi(:,outID), Interaction, InteractionFraction, ...
        'VariableNames', {'Parameter','Si','STi','Interaction','InteractionFraction'});

    T = sortrows(T,'STi','descend');

    writetable(T,filename, 'Sheet',outputNames{outID});

end

plot_sobol(Si,STi,parNames,outputNames,N);

%% Export summary
threshold = 0.1; % count outputs where the parameter explains at least 10% of the variance
peakIdx = 1:8;
ttpIdx  = 9:15;
filename = sprintf('output/Sobol_global_summary_N%d.xlsx',N);

%% ALL OUTPUTS
MeanSTi = mean(STi,2);
MaxSTi  = max(STi,[],2);
Count   = sum(STi > threshold,2);

RelevantOutputs = strings(length(parNames),1);

for p = 1:length(parNames)
    idx = find(STi(p,:) > threshold);
    if isempty(idx)
        RelevantOutputs(p) = "";
    else
        RelevantOutputs(p) = strjoin(outputNames(idx), ', ');
    end
end

SummaryAll = table( ...
    parNames(:), ...
    MeanSTi, ...
    MaxSTi, ...
    Count, ...
    RelevantOutputs, ...
    'VariableNames',{ ...
    'Parameter', ...
    'MeanSTi', ...
    'MaxSTi', ...
    'NumOutputs', ...
    'RelevantOutputs'});

SummaryAll = sortrows(SummaryAll,'MeanSTi','descend');
writetable(SummaryAll,filename,'Sheet','AllOutputs');

%% PEAK OUTPUTS
MeanSTi = mean(STi(:,peakIdx),2);
MaxSTi  = max(STi(:,peakIdx),[],2);
Count   = sum(STi(:,peakIdx) > threshold,2);

RelevantOutputs = strings(length(parNames),1);

for p = 1:length(parNames)
    idx = peakIdx(STi(p,peakIdx) > threshold);
    if isempty(idx)
        RelevantOutputs(p) = "";
    else
        RelevantOutputs(p) = strjoin(outputNames(idx), ', ');
    end
end

SummaryPeak = table( ...
    parNames(:), ...
    MeanSTi, ...
    MaxSTi, ...
    Count, ...
    RelevantOutputs, ...
    'VariableNames',{ ...
    'Parameter', ...
    'MeanSTi', ...
    'MaxSTi', ...
    'NumOutputs', ...
    'RelevantOutputs'});

SummaryPeak = sortrows(SummaryPeak,'MeanSTi','descend');
writetable(SummaryPeak,filename,'Sheet','Peaks');

%% TIME-TO-PEAK OUTPUTS
MeanSTi = mean(STi(:,ttpIdx),2);
MaxSTi  = max(STi(:,ttpIdx),[],2);
Count   = sum(STi(:,ttpIdx) > threshold,2);

RelevantOutputs = strings(length(parNames),1);

for p = 1:length(parNames)
    idx = ttpIdx(STi(p,ttpIdx) > threshold);
    if isempty(idx)
        RelevantOutputs(p) = "";
    else
        RelevantOutputs(p) = strjoin(outputNames(idx), ', ');
    end
end

SummaryTTP = table( ...
    parNames(:), ...
    MeanSTi, ...
    MaxSTi, ...
    Count, ...
    RelevantOutputs, ...
    'VariableNames',{ ...
    'Parameter', ...
    'MeanSTi', ...
    'MaxSTi', ...
    'NumOutputs', ...
    'RelevantOutputs'});

SummaryTTP = sortrows(SummaryTTP,'MeanSTi','descend');
writetable(SummaryTTP,filename,'Sheet','TimeToPeak');

% Finish
elapsedTime = toc;
fprintf('\nFinished %d simulations.\n',nSim);
fprintf('Elapsed time: %.2f seconds\n',elapsedTime);
fprintf('Average time per simulation: %.4f seconds\n',elapsedTime/nSim);