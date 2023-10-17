%% This is used for calculate average spectrum and summary statistics

% 2023/10/07 By Zijiang Yang
% 2023/10/09 Adding CI calculations

%% Prepare workspace
clc, clear, close all

tic

%% Input data
% Relative path to the 'Input data' folder
filename     = 'A02_std spectra.csv';

relativePath = '..\A01_Data\';

% Construct the full path for the data files
Namesmp = [relativePath, filename];

% Open the file for reading
fid = fopen(Namesmp, 'r');

% Read the first line and process "Shape"
line = fgetl(fid);
tokens = split(line, ',');
Shape = string(tokens(2:end));

% Read the second line and process "Color"
line = fgetl(fid);
tokens = split(line, ',');
Color = string(tokens(2:end));

% Close the file
fclose(fid);

% Read the remaining data into a table starting from the 3rd row
opts = detectImportOptions(Namesmp, 'NumHeaderLines', 2);
dataTable  = readtable(Namesmp, opts);

clear fid line tokens

%% Parameters
para.wd = 25;       % DSW method window size. If para.wd = 0, then not used. Window size; may be 10:5:200 for optimization; window size and step size for baseline Lenz et al., 2015 para.ws = 128 para.ss = 128; para.wd and 5*para.wd are the two window sizes for two potential baselines. The final estimated baseline is based on these two potetial baselines
para.hq = 0.4;      % Default 0.4/ flowcell with ethanol 0.6/HQI ranges from -1 to 1; may be 0.0:0.05:0.8 
para.dn = 0;        % direvative order: para.dn = 0 no; para.dn = 1 1st order; para.dn = 2 2nd order; 1st order seemed to be the best (Renner et al., 2019)
para.di = 0.00001;  % co-polymer shreshold, default 0.1

para.dv  = 0.5;         % default 0.5 alignment interval; accuracy of spectrum 1D-interpolation
para.sp  = 11;          % smoothing range, Nakano et al., 2021 para.sp = 11
para.qt  = 0;           % Quantile for localized baseline; Lenz et al., 2015 para.qt = 0
para.pr  = 5000;         % moving windowsize for pruning, 25 is an empirical value
para.rg  = [400 4000];  % range of measured spectrum, Nakano et al., 2021 para.sp = [400 4000]
para.ct  = [0 1]; % For ATR [650 700; 2250 2450; 000 500; 3500 8000], Nakano et al., 2021 
para.sm  = 1;           % para.sm = 0 for moving average (Nakano et al., 2021); para.sm = 1 for Savitzky-Golay (Lenz et al., 2015)
para.od  = 2;           % polynomial order for Savitzky-Golay

para.X1 = [1650 1850];  % X1, X2 = if a wavenumber1 and wavenumber2, then using area method; if wavenumber and 5, then using peak height method
para.X2 = [1420 1500];  % e.g.., X1 = [1650 1850] or X1 = [1650 5]
para.method = 'G';      % method: baseline method. if = 'G', then global baseline correction, if = 'L', then local baseline correction (Walfridson and Kuttainen, 2022)

%% Spectral processing
% Spectral preprocessing (I) - baseline correction
[SNR,pre1XY,~,~]  = Function_DSW_f04(table2array(dataTable),para.wd,0);

% Spectral proprocessing (II) - pruning
[pre2XY0,pre2XY1,pre2XY2] = Function_Pruning_f02(pre1XY,para.dv,para.sp,0,para.qt,para.pr,para.pr,para.rg,para.ct,para.sm,para.od);

% Spectral proprocessing (II) - normalization
min_val_XY        = min(pre2XY0(:, 2));              % 1. Find the minimum value of the second column
pre2XY0(:, 2:end) = pre2XY0(:, 2:end) - min_val_XY;  % 2. Subtract this minimum value from all columns from the 2nd column onward
max_val_XY        = max(pre2XY0(:, 2:end));          % 3. Find the maximum value of the updated second column
pre2XY0(:, 2:end) = pre2XY0(:, 2:end) ./ max_val_XY; % 4. Divide all columns (from the 2nd column onward) by this maximum value

% Organize results
X = pre2XY0(:,1);
Y = pre2XY0(:,2:end);

%% Color-shape counting
% 1. Get unique colors and shapes
uniqueColors = unique(Color);
uniqueShapes = unique(Shape);

numUniqueColors = numel(uniqueColors);
numUniqueShapes = numel(uniqueShapes);

% 2. Initialize table with zeros
tableCounts = zeros(numUniqueColors, numUniqueShapes);

% 3. Populate the table
for i = 1:numUniqueColors
    for j = 1:numUniqueShapes
        tableCounts(i,j) = sum(Color == uniqueColors(i) & Shape == uniqueShapes(j));
    end
end

% 4. Compute the SUM values
sumColors = sum(tableCounts, 2); % sum across columns
sumShapes = sum(tableCounts, 1); % sum across rows

% 5. Add SUM values to the table
tableCounts = [tableCounts, sumColors];
tableCounts = [tableCounts; [sumShapes, sum(sumColors)]];

% 6. Create a better reporting format using the MATLAB table datatype
T.count = array2table(tableCounts, 'RowNames', [cellstr(uniqueColors); 'SUM'], 'VariableNames', [cellstr(uniqueShapes); 'SUM']);

% Display the table
disp(T.count);

%% Averaging spectra I - seperation
% Create matrices based on unique Color
for i = 1:numUniqueColors
    currentColor = uniqueColors{i}; % Get the color name
    indices = strcmp(Color, currentColor); % Find the indices matching this color
    Ycolor.(currentColor) = Y(:, indices); % Store the matching columns in the Ycolor structure
end

% Create matrices based on unique Shape
for j = 1:numUniqueShapes
    currentShape = uniqueShapes{j}; % Get the shape name
    indices = strcmp(Shape, currentShape); % Find the indices matching this shape
    Yshape.(currentShape) = Y(:, indices); % Store the matching columns in the Yshape structure
end

% Create matrices based on the color-shape combination
for i = 1:numUniqueColors
    for j = 1:numUniqueShapes
        currentColor = uniqueColors{i}; % Get the color name
        currentShape = uniqueShapes{j}; % Get the shape name
        combinationName = [currentColor, '_', currentShape]; % Create the combination name
        indices = strcmp(Color, currentColor) & strcmp(Shape, currentShape); % Find the indices matching this combination
        Ycosh.(combinationName) = Y(:, indices); % Store the matching columns in the Ycosh structure
    end
end

%% Averaging spectra II - mean
% Calculate the mean for each matrix in the Ycolor structure
fieldsColor = fieldnames(Ycolor);
for i = 1:length(fieldsColor)
    fieldName = fieldsColor{i};
    Yavg.mean.color.(fieldName) = mean(Ycolor.(fieldName), 2);
end

% Calculate the mean for each matrix in the Yshape structure
fieldsShape = fieldnames(Yshape);
for i = 1:length(fieldsShape)
    fieldName = fieldsShape{i};
    Yavg.mean.shape.(fieldName) = mean(Yshape.(fieldName), 2);
end

% Calculate the mean for each matrix in the Ycosh structure:
fieldsCosh = fieldnames(Ycosh);
for i = 1:length(fieldsCosh)
    fieldName = fieldsCosh{i};
    Yavg.mean.cosh.(fieldName) = mean(Ycosh.(fieldName), 2);
end

% Finally, you can generate the desired table:
tableColorMeans = struct2table(Yavg.mean.color);
tableShapeMeans = struct2table(Yavg.mean.shape);
tableCoshMeans = struct2table(Yavg.mean.cosh);
YavgTableMean = [tableColorMeans, tableShapeMeans, tableCoshMeans];

%% Averaging spectra III - 5% quantile
fieldsColor = fieldnames(Ycolor);
for i = 1:length(fieldsColor)
    fieldName = fieldsColor{i};
    Yavg.Q05.color.(fieldName) = quantile(Ycolor.(fieldName), 0.05, 2);
end

fieldsShape = fieldnames(Yshape);
for i = 1:length(fieldsShape)
    fieldName = fieldsShape{i};
    Yavg.Q05.shape.(fieldName) = quantile(Yshape.(fieldName), 0.05, 2);
end

fieldsCosh = fieldnames(Ycosh);
for i = 1:length(fieldsCosh)
    fieldName = fieldsCosh{i};
    Yavg.Q05.cosh.(fieldName) = quantile(Ycosh.(fieldName), 0.05, 2);
end

% Convert the Q05 matrices from each structure into tables
tableColorQ05 = struct2table(Yavg.Q05.color);
tableShapeQ05 = struct2table(Yavg.Q05.shape);
tableCoshQ05 = struct2table(Yavg.Q05.cosh);

% Combine these tables
YavgTableQ05 = [tableColorQ05, tableShapeQ05, tableCoshQ05];

%% Averaging spectra IV - 95% quantile
fieldsColor = fieldnames(Ycolor);
for i = 1:length(fieldsColor)
    fieldName = fieldsColor{i};
    Yavg.Q95.color.(fieldName) = quantile(Ycolor.(fieldName), 0.95, 2);
end
fieldsShape = fieldnames(Yshape);
for i = 1:length(fieldsShape)
    fieldName = fieldsShape{i};
    Yavg.Q95.shape.(fieldName) = quantile(Yshape.(fieldName), 0.95, 2);
end
fieldsCosh = fieldnames(Ycosh);
for i = 1:length(fieldsCosh)
    fieldName = fieldsCosh{i};
    Yavg.Q95.cosh.(fieldName) = quantile(Ycosh.(fieldName), 0.95, 2);
end
% Convert the Q95 matrices from each structure into tables
tableColorQ95 = struct2table(Yavg.Q95.color);
tableShapeQ95 = struct2table(Yavg.Q95.shape);
tableCoshQ95 = struct2table(Yavg.Q95.cosh);

% Combine these tables
YavgTableQ95 = [tableColorQ95, tableShapeQ95, tableCoshQ95];

%% CI calculation
% Calculate the CI for matrices in Ycolor
fieldsColor = fieldnames(Ycolor);
for i = 1:length(fieldsColor)
    fieldName = fieldsColor{i};
    YY = Ycolor.(fieldName);
    for j = 1:size(YY, 2)
        CI.Color.(fieldName)(:, j) = Function_CI_calculator_f05([X YY(:, j)], para.X1, para.X2, para.method);
    end
end

% Calculate the CI for matrices in Yshape
fieldsShape = fieldnames(Yshape);
for i = 1:length(fieldsShape)
    fieldName = fieldsShape{i};
    YY = Yshape.(fieldName);
    for j = 1:size(YY, 2)
        CI.Shape.(fieldName)(:, j) = Function_CI_calculator_f05([X YY(:, j)], para.X1, para.X2, para.method);
    end
end

% Calculate the CI for matrices in Ycosh
fieldsCosh = fieldnames(Ycosh);
for i = 1:length(fieldsCosh)
    fieldName = fieldsCosh{i};
    YY = Ycosh.(fieldName);
    for j = 1:size(YY, 2)
        CI.ShCo.(fieldName)(:, j) = Function_CI_calculator_f05([X YY(:, j)], para.X1, para.X2, para.method);
    end
end


%% CI comparison
categories = {'Shape', 'Color', 'ShCo'};
pTables = struct();

for cat = categories
    category = cat{1};

    % Log10 Transformation
    fields = fieldnames(CI.(category));
    for i = 1:length(fields)
        fieldName = fields{i};
        CI.([category 'Log']).(fieldName) = log10(CI.(category).(fieldName));
    end
    
    % Welch Tests
    fieldsLog = fieldnames(CI.([category 'Log']));
    numFields = length(fieldsLog);
    
    pValuesMatrix = ones(numFields); % Initializing to ones
    for i = 1:numFields
        for j = 1:numFields
            if i ~= j % don't compare data with itself
                [~, p] = ttest2(CI.([category 'Log']).(fieldsLog{i}), CI.([category 'Log']).(fieldsLog{j}), 'Vartype','unequal');
                pValuesMatrix(i,j) = p;
            end
        end
    end
    
    % Convert P-values
    pValuesStrMatrix = cell(size(pValuesMatrix));
    pValuesStrMatrix(pValuesMatrix < 0.05) = {'diff'};
    pValuesStrMatrix(pValuesMatrix >= 0.05) = {'same'};
    
    % Convert to table
    pTable = array2table(pValuesStrMatrix, 'RowNames', fieldsLog, 'VariableNames', fieldsLog);
    
    % Save to the structured output
    pTables.(category) = pTable;
end




% Summarize CI values:
% For Shape
fieldsShape = fieldnames(CI.Shape);
maxLengthShape = max(structfun(@(x) size(x,1), CI.Shape));  % This line was missing
CI_Summary_Shape = NaN(maxLengthShape, length(fieldsShape));

for i = 1:length(fieldsShape)
    currentData = CI.Shape.(fieldsShape{i})(:);
    currentLength = length(currentData);
    CI_Summary_Shape(1:currentLength, i) = currentData;
end

% For Color
fieldsColor = fieldnames(CI.Color);
maxLengthColor = max(structfun(@(x) size(x,1), CI.Color));
CI_Summary_Color = NaN(maxLengthColor, length(fieldsColor));

for i = 1:length(fieldsColor)
    currentData = CI.Color.(fieldsColor{i})(:);
    currentLength = length(currentData);
    CI_Summary_Color(1:currentLength, i) = currentData;
end

% For ShCo
fieldsShCo = fieldnames(CI.ShCo);
maxLengthShCo = max(structfun(@(x) size(x,1), CI.ShCo));
CI_Summary_ShCo = NaN(maxLengthShCo, length(fieldsShCo));

for i = 1:length(fieldsShCo)
    currentData = CI.ShCo.(fieldsShCo{i})(:);
    currentLength = length(currentData);
    CI_Summary_ShCo(1:currentLength, i) = currentData;
end

% Ensure all matrices have the same number of rows
maxLengthShape = height(CI_Summary_Shape);
maxLengthColor = height(CI_Summary_Color);
maxLengthShCo = height(CI_Summary_ShCo);

maxLength = max([maxLengthShape, maxLengthColor, maxLengthShCo]);

CI_Summary_Shape(end+1:maxLength, :) = NaN;
CI_Summary_Color(end+1:maxLength, :) = NaN;
CI_Summary_ShCo(end+1:maxLength, :) = NaN;

CI_Summary_Shape = [CI_Summary_Shape; NaN(maxLength - size(CI_Summary_Shape, 1), size(CI_Summary_Shape, 2))];
CI_Summary_Color = [CI_Summary_Color; NaN(maxLength - size(CI_Summary_Color, 1), size(CI_Summary_Color, 2))];
CI_Summary_ShCo = [CI_Summary_ShCo; NaN(maxLength - size(CI_Summary_ShCo, 1), size(CI_Summary_ShCo, 2))];

% Concatenate
CI_Summary = [CI_Summary_Shape, CI_Summary_Color, CI_Summary_ShCo];


% Convert to a table
CI.Summary = array2table(CI_Summary, 'VariableNames', [fieldsShape; fieldsColor; fieldsShCo]);

% Iterate over all variables (columns) in the table
for varName = CI.Summary.Properties.VariableNames
    currentVar = varName{1}; % Because VariableNames is a cell array
    CI.Summary.(currentVar)(CI.Summary.(currentVar) == 0) = NaN;
end



%% The final outputs are:
Results.Tcount = T.count;
Results.X      = array2table(X, 'VariableNames', {'X'});
Results.YavgTableMean = YavgTableMean;
Results.YavgTableQ05 = YavgTableQ05;
Results.YavgTableQ95 = YavgTableQ95;

Results.CIsummary = CI.Summary;
Results.p_Shape = pTables.Shape;
Results.p_Color = pTables.Color;
Results.p_ShCo = pTables.ShCo;


% Extract name from your given filename
[~, namePart, ~] = fileparts(filename);

% Create a new filename with the desired format
currentDateTime = datetime('now', 'Format', 'yyyyMMddHHmm');
newFilename = sprintf('%s_%s.xlsx', namePart, currentDateTime);

% Iterate over the fields in Results and write each table to the Excel file
fields = fieldnames(Results);
for i = 1:length(fields)
    tbl = Results.(fields{i});
    writetable(tbl, newFilename, 'Sheet', fields{i});
end

% Create two cell arrays for parameter names and values
paramNames = fieldnames(para);
paramValues = struct2cell(para);

% Convert these cell arrays to a table
paramTable = table(paramNames, paramValues);
paramTable.Properties.VariableNames = {'Parameter', 'Value'};

% Write the parameter table to the Excel file
writetable(paramTable, newFilename, 'Sheet', 'parameters');

% Get all table names in the Results structure
tableNames = fieldnames(Results);

for i = 1:length(tableNames)
    currentTable = Results.(tableNames{i});
    
    % If the table has row names, convert them to a new column
    if ~isempty(currentTable.Properties.RowNames)
        rowNames = currentTable.Properties.RowNames;
        currentTable = addvars(currentTable, rowNames, 'Before', 1, 'NewVariableNames', 'Items');
    end
    
    % Write to the Excel file
    writetable(currentTable, newFilename, 'Sheet', tableNames{i});
end


%% Ending
% Convert the elapsed time to minutes and seconds
elapsedTime = toc;  % This gets the elapsed time in seconds
minutes = floor(elapsedTime/60);
seconds = rem(elapsedTime, 60);
message = sprintf('ご主人様、計算は終りましたにゃ。計算時間は%d分%.1f秒でしたにゃ', minutes, seconds);

% Display the message box
msgbox(message, 'Calculation Completed');

toc



