function data = getData(folder, theFiles, condition)
%  ------------------------------------------------------------------------
%  Author: Siddhartha Dhiman
%  E-Mail: sdhiman@buffalo.edu
%  Created: 07/24/2017 using Matlab R2016b
%  Credits to David Legland for Euclidean distances and MST, and matGeom
%  toolbox for graphing
%  -------------------------------------------------------------------------
%  Function to generate minimum distances and coordinates into the matrix
%  structure 'data'. This matrix contains the following information:
%       FileName:   Name of files processed
%       ROI:        ROI index number
%       MinDist:    Vector containing all Euclidean minimum distances
%                   computed from MST
%       Mean:       Average of MinDist
%       StDev:      Standard deviation of MinDist
%   -----------------------------------------------------------------------
%   Input Arguments:
%       folder:     Path to folder containing images
%       theFiles:   Structural matrix with directory information. Can be
%                   constructed using 'dir(filePattern)', where
%                   'filePattern' is a file filtering criteria.
%       condition:  1 for Control
%                   2 for Schizophrenia
%   -----------------------------------------------------------------------
%   Output:
%       getData: Structural matrix containing information useful for
%                generating biological networks.
%   -----------------------------------------------------------------------
%   Feel free to drop me an email on questions and concerns
%   -----------------------------------------------------------------------

%% Variables for Cell Detection
%   These variables can be changed to fine tune the image output
%   -----------------------------------------------------------------------
R_SE1 = 15;             % Top-Hat structuring element radius (counting)
R_SE2 = 5;              % Morphological opening structuring element radius
sq = 5;                 % Size of 'sq x sq' square matrix for marker cleaning
open_AC = 20;           % Area for noise removal using 'bwareaopen'
H = 5;                  % H-maxima transform
minCellSize = 10;       % Minimum cell size in Watershed transform
sz = 25;                % Scatter plot point size

%% Structuring Elements for Cell Detection
SE1 = strel('disk', R_SE1);     % For Top-Hat (Counting)
SE2 = strel('disk', R_SE2);     % For morphological operations
SE3 = strel(ones(sq,sq));       % For marker cleaning

%% Initialize Data Matrix
field1 = 'FileName';
field2 = 'ROI';
field3 = 'CellCount';
field4 = 'Nodes';
field5 = 'MinDist';
field6 = 'Mean';
field7 = 'StDev';
value1 = 0;
value2 = 0;
value3 = 0;
value4 = 0;
value5 = 0;
value6 = 0;
value7 = 0;
data = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
    field5,value5,field6,value6,field7,value7);

%% Setting up Loop Counter
cntA = 0;
cntB = 0;

%% Main Loop
for a = 1:length(theFiles)
    cntA = cntA + 1;
    baseFileName = theFiles(a).name;
    fullFileName = fullfile(folder, baseFileName);
    fprintf(1, 'Now evaluating %s\n', baseFileName);
    imageArray = imread(fullFileName);
    [pathstr fileName ext] = fileparts(fullFileName);
    
    %   Extract red channel
    IR = imageArray(:,:,1);
    
    %% Top-Hat Operation
    tophatIR = imtophat(IR, SE1);
    eqIR = adapthisteq(tophatIR);
    
    %% Morphological Processing
    IRo = imopen(eqIR, SE2);
    erodeIR = imerode(tophatIR, SE2);
    IRobr = imreconstruct(erodeIR, tophatIR);
    IRoc = imclose(IRo, SE2);
    IRobrd = imdilate(IRobr, SE2);
    IRobrcbr = imreconstruct(imcomplement(IRobrd), imcomplement(IRobr));
    IRobrcbr = imcomplement(IRobrcbr);
    
    %% Regional Maximas
    fgm = imextendedmax(IRobrcbr, H);
    I2 = imageArray;
    I2(fgm) = 255;
    fgm2 = imclose(fgm, SE3);
    fgm3 = imerode(fgm2, SE3);
    fgm4 = bwareaopen(fgm3, open_AC);
    I3 = imageArray;
    I3(fgm4) = 255;
    bw = imbinarize(IRobrcbr);
    bw4_perim = bwperim(bw);
    
    %% Watershed Transform for Cell Detection
    eqIR_c = imcomplement(eqIR);
    I_mod = imimposemin(eqIR_c, ~bw | fgm4);
    L = watershed(I_mod);
    
    binWatershed = L > 1;
    regsIR = regionprops(binWatershed, 'Area', 'Centroid', 'PixelIdxList');
    
    %% ROI Selection and Masking
    R_circ = 5;         % Radius of re-plotted circle
    imageSizeX = size(IR,2);
    imageSizeY = size(IR,1);
    IRcond = zeros(imageSizeY, imageSizeX);
    [columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
    
    ROI = ROIdetection(imageArray);
    sizeROI = size(ROI,3);
    for b = 1:sizeROI
        cntB = cntB + 1;
        data(cntB).FileName = baseFileName;
        data(cntB).ROI = b;
        ROImask = ROI(:,:,b);
        
        for k = 1:numel(regsIR)     % Replots all detected cells into circles
            circles = (rowsInImage - regsIR(k).Centroid(2)).^2 ...
                + (columnsInImage - regsIR(k).Centroid(1)).^2 <= R_circ.^2;
            IRcond = IRcond + circles;
        end
        
        cellsROI = imbinarize(IRcond .* ROImask);
        
        %% Statistics from ROI
        regsCells = regionprops(cellsROI, 'Centroid', 'Extrema');
        data(cntB).CellCount = numel(regsCells);
        
        %   Position is an NxM matrix containing node coordinates where N is the
        %   number of cells detected in ROI and M is the number of dimensions
        %   In this case, M = 2 so this is a two-dimensional problem
        
        for i = 1:numel(regsCells)
            Position(i,1) = regsCells(i).Centroid(1);   % Create X-coordinate row vector
            Position(i,2) = regsCells(i).Centroid(2);   % Create Y-coordinate row vector
        end
        %   Saves node coordinates into data matrix
        data(cntB).Nodes = Position;
        
        %% Compute Euclidean Distances
        D   = size(Position, 2);
        Df  = factorial(D);
        
        %   Find all possible vertices in unit triangulation
        %   Variable 'subs' shows all possible combinations
        subs = zeros(Df, 2);
        k = 1;
        for i = 1:D
            for j = i+1:D+1
                subs(k, 1) = i;
                subs(k, 2) = j;
                k = k + 1;
            end
        end
        
        %   Compute Delaunay Triangulation
        TRI = delaunayn(Position);
        tsize = size(TRI, 1);
        
        %   Compute all possible edges
        edges = zeros(tsize*Df, 2);
        for t = 1:tsize
            for i = 1:Df
                edges((t-1)*Df+i, 1) = TRI(t, subs(i, 1));
                edges((t-1)*Df+i, 2) = TRI(t, subs(i, 2));
            end
        end
        
        %   Simplify Edges
        edges = unique(sort(edges, 2), 'rows');
        
        %   Compute Euclidean distance of each edge
        dist = zeros(size(edges, 1), 1);
        for i = 1:size(edges,1)
            dist(i) = distancePoints(Position(edges(i,1), :), Position(edges(i,2), :));
        end
        
        %% Form MST using Prim's Algorithm
        %   Isolate vertices index
        nodes = unique(edges(:));
        N = length(nodes);
        
        %   Initialize memory
        nodes2 = zeros(N, 1);
        EDGES = zeros(N-1, 2);
        eucD = zeros(N-1, 1);
        
        %   Initialize with a first node
        nodes2(1) = nodes(1);
        nodes = nodes(2:end);
        
        %   Iterate on edges
        for i = 1:N-1
            %   Find all edges from nodes2 to nodes
            ind = unique(find(...
                (ismember(edges(:,1), nodes2(1:i)) & ismember(edges(:,2), nodes)) | ...
                (ismember(edges(:,1), nodes) & ismember(edges(:,2), nodes2(1:i))) ));
            
            %   Choose edge with lowest value
            [tmp, ind2] = min(dist(ind)); %#ok<ASGLU>
            ind = ind(ind2(1));
            eucD(i) = dist(ind);
            
            %   Index of other vertex
            edge = edges(ind, :);
            neigh = edge(~ismember(edge, nodes2));
            
            %   Add to list of nodes and list of edges
            nodes2(i+1) = neigh;
            EDGES(i,:) = edge;
            
            %   Remove current node from list of old nodes
            nodes = nodes(~ismember(nodes, neigh));
        end
        
        data(cntB).MinDist = eucD;
        
        %% Generate MSTs
        
        gcf = figure('vis', 'off');
        drawGraph(Position, EDGES)
        xlabel('Distance (Pixel)');
        ylabel('Distance (Pixel)');
        if condition == 1
            saveFolder = fullfile(pwd, 'Output', 'MSTs','Control', fileName);
            title('Control Ki-16 MST');
        elseif condition == 2
            saveFolder = fullfile(pwd, 'Output','MSTs', 'Schizophrenia', fileName);
            title('Schizophrenia Ki-16 MST');
        else
            msg = 'Condition not specified or invalid when calling getData';
            error(msg)
        end
        
        mkdir(saveFolder);
        print(gcf, fullfile(saveFolder, sprintf('ROI_%d',b)),'-dpng');

        %% Generate Supplemental Material
        hbins = 6;
        green = cat(3, zeros(size(IR)),ones(size(IR)), zeros(size(IR)));
        gch = figure('vis', 'off');
        subplot(2,1,1);
        imshow(imageArray), hold on;
        h = imshow(green); 
        hold off;
        set(h, 'AlphaData', cellsROI);
        axis image;
        title('Detected Cells Imposed on Original Image');
        
        subplot(2,1,2);
        h = histogram(eucD, hbins);
        title('Histogram of Euclidean Distances');
        
        if condition == 1
            saveFolder = fullfile(pwd, 'Output', 'Supplemental','Control', fileName);
        elseif condition == 2
            saveFolder = fullfile(pwd, 'Output','Supplemental', 'Schizophrenia', fileName);
        else
            msg = 'Condition not specified or invalid when calling getData';
            error(msg)
        end
        mkdir(saveFolder);
        print(gch, fullfile(saveFolder, sprintf('ROI_%d',b)),'-dpng');
        
        %% Calculate Mean and Deviation
        data(cntB).Mean = mean(eucD);
        data(cntB).StDev = std(eucD);
        
        %% Delete Loop Variables for Repeating
        clearvars regsCells Position eucD
        
    end
end

data = data;
clearvars variables -except data