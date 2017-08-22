%%  MATLAB Script for Counting Cells and Forming MST
%  ------------------------------------------------
%  Author: Siddhartha Dhiman
%  E-Mail: sdhiman@buffalo.edu
%  Created: 07/24/2017 using Matlab R2016b
%  Credits to David Legland for Euclidean distances and MST, and matGeom
%  toolbox for graphing

%% Clearing Workspace
clc;
clear all;
close all;
warning off;

%% Select Folder
folderCtrl = uigetdir('', 'Select Control Image Directory');
if ~isdir(folderCtrl)
    errorMessage = sprintf('Error: The following directory does not exist:\n%s', folderCtrl);
    uiwait(warndlg(errorMessage));
    return;
end

filePatternCtrl = fullfile(folderCtrl, '*.tif');
theFilesCtrl = dir(filePatternCtrl);
[upperPathCtrl, deepestFolderCtrl, ~] = fileparts(folderCtrl);

folderSch = uigetdir('', 'Select Schizophrenia Image Directory');
if ~isdir(folderSch)
    errorMessage = sprintf('Error: The following directory does not exist:\n%s', folderSch);
    uiwait(warndlg(errorMessage));
    return;
end

filePatternSch = fullfile(folderSch, '*.tif');
theFilesSch = dir(filePatternSch);
[upperPathCtrlSch, deepestFolderCtrlSch, ~] = fileparts(folderSch);

%% Generate Data Matrix
dataCtrl = getData(folderCtrl, theFilesCtrl, 1);
fprintf(1, '\n');
fprintf(1, ' Data matrix dataCtrl (Control) created \n');
fprintf(1, '\n');
fprintf(1, '\n');

dataSch = getData(folderSch, theFilesSch, 2);
fprintf(1, '\n');
fprintf(1, 'Data matrix dataSch (Schizophrenia) created \n');
fprintf(1, '\n');
fprintf(1, '\n');

%% Vectorize Distances
%   Control Group
cVector = vertcat(dataCtrl(:).MinDist);
    
%   Schizophrenia Group
sVector = vertcat(dataSch(:).MinDist);

%% Calculate Global Mean and St. Dev
muC = mean([dataCtrl.Mean]);
muS = mean([dataSch.Mean]);
stdevC = std([dataCtrl.Mean]);
stdevS = std([dataSch.Mean]);

nbins = 50;
gcf = figure;
hC = histogram(cVector,nbins);
hC.FaceColor = 'r';
hC.EdgeAlpha = 0.5;
hC.EdgeColor = 'none';
hold on;
hS = histogram(sVector, nbins);
hS.FaceColor = 'b';
hS.EdgeAlpha = 0.5;
hS.EdgeColor = 'none';
hold off;
legend('Control','Schizophrenia')
title('Histograms of Global Control and Schizophrenia Groups')
saveFolder = fullfile(pwd, 'Output','Global');
mkdir(saveFolder);
print(gcf, fullfile(saveFolder, sprintf('Global Histogram')),'-dpng');

%% Export to CSV
rm = {'Nodes', 'MinDist'};
dataC = struct2table(rmfield(dataCtrl, rm));
dataS = struct2table(rmfield(dataSch, rm));

writetable(dataC, fullfile(saveFolder, 'Control.csv'));
writetable(dataS, fullfile(saveFolder, 'Schizophrenia.csv'));

%% Binning Per Group
T = [cVector',sVector'];
a = min(T);
a = round(a,-1);
b = max(T);
b = round(b,-1);
space = 51;
%   Create binning edges
d = linspace(a,b,space);

for i = 1:(length(d)-1)
    clear s t
    s  = find(cVector >= d(i) & cVector < d(i+1));
    t = find(sVector >= d(i) & sVector < d(i+1));
    vc = cVector(s);
    vs = sVector(t);
    binMat(i).Control = vc;
    binMat(i).Schizophrenia = vs;
    
%     if any(binMat(i).Control) == 0
%         binMat(i).Control = NaN;
%     else
%         binMat(i).Control = binMat(i).Control;
%     end
%     
%     if any(binMat(i).Schizophrenia) == 0
%         binMat(i).Schizophrenia = NaN;
%     else
%         binMat(i).Schizophrenia = binMat(i).Schizophrenia;
%     end
end

writetable(struct2table(binMat),fullfile(saveFolder,'Histogram Array.xlsx'));

%% Frequency Table
for i = 1:length(binMat)
    BinID{i,1} = ['Bin ' num2str(i)];
    Control(i,1) = length(binMat(i).Control);
    Schizophrenia(i,1) = length(binMat(i).Schizophrenia);
end

tableFreq = table(Control, Schizophrenia, 'RowNames', BinID);
writetable(tableFreq, fullfile(saveFolder, 'Frequency Table.csv'),...
    'WriteRowNames',true);


%% Binning Per ROI
%   Initiliaze Structuring Arrays
for i = 1:length(dataCtrl)
    headerROIc(i) = {sprintf('ROI_%d',i)};
end
headerROIc = char(headerROIc);
binROIc.(headerROIc) = 0;
    
    
for i = 1:length(dataCtrl)
    binROIc = 0;
end

for i = 1:(length(d)-1)
    for k = 1:length(dataCtrl)
        clear s vc
        s = find(dataCtrl(k).MinDist >= d(i) & dataCtrl(k).MinDist ...
            < d(i+1));
        
        vc = dataCtrl(k).MinDist(s);
        
        binROIc(k,i) = length(vc);
        
    end
    for l = 1:length(dataSch)
        clear t vs
        t = find(dataSch(l).MinDist >= d(i) & dataSch(l).MinDist ...
            < d(i+1));
        vs = dataSch(l).MinDist(t);
        binROIs(l,i) = length(vs);
    end
end

    
    

%% Chi-Sqared Testing
%   Initialize Chi Matrix
nC = sum(Control);
nS = sum(Schizophrenia);
for i = 1:length(binMat);
    clear x1 x2 tbl chi2stat pval
    x1 = [repmat('a',nC,1); repmat('b',nS,1)];
    x2 = [repmat(1,Control(i),1); repmat(2,nC-Control(i),1); ... 
        repmat(1,Schizophrenia(i),1); repmat(2,nS-Schizophrenia(i),1)];
    [tbl, chi2stat, pval] = crosstab(x1,x2);
    chi(i,1) = chi2stat;
    pvalue(i,1) = pval;
end
tableChi = table(chi, pvalue, 'RowNames', BinID);
writetable(tableChi, fullfile(saveFolder, 'Chi Test.xlsx'),...
    'WriteRowNames',true);

    