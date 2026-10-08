%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 17/09/2026                  %
% ----------------------------------------------------------- %
%  Post-processing of the results after running the file      %
%  Sobol_convergence.m. Files needed:                         %
%    * plot_sobol.m                                           %
%    * Sobol_convergence.mat                                  %
% Additional main codes:                                      %
%    * Sobol_convergence.m: check convergence to choose N     %
%    * Sobol_main.m: runs global sensitivity analaysis        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('output-convergence/Sobol_convergence.mat')

% Choose among Nlist = [64 128 256 512 1024 2048 4096 8192]
N = 8192;

nID = find([Results.N] == N);

Si  = Results(nID).Si;
STi = Results(nID).STi;

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
ttpIdx  = 9:16;
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