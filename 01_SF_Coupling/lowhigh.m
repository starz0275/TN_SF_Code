% =========================================================================
% Structural-Functional (S-F) Coupling Calculation using Graph Harmonic Model
% =========================================================================
clear all; clc;

% -------------------------------------------------------------------------
% 1. Set Directory Paths (Please update these relative paths to your local data)
% -------------------------------------------------------------------------
FCdirpath = './Data_Sample/Fmri_fc';      % Path to functional connectivity matrices
SCdirpath = './Data_Sample/DTI_data_sym'; % Path to structural connectivity matrices

FCfiles = dir(FCdirpath);                 % List all files in the FC directory
SCfiles = dir(SCdirpath);                 % List all files in the SC directory

% Initialize variables
allrsq = zeros(400, length(FCfiles)-2);
allxx = zeros(400, length(FCfiles)-2);
numflag = 0;

% -------------------------------------------------------------------------
% 2. Main Loop: Iterate through subjects and match multimodal data
% -------------------------------------------------------------------------
for ii = 1:length(FCfiles)
    flagsc = 0;
    flagfc = 0;
    flagmn = 1; % Status flag
    
    % Skip directories (. and ..)
    if FCfiles(ii).isdir
        continue;
    end
    
    fc_file_name = FCfiles(ii).name;
    
    % Extract FC subject ID
    tokens = regexp(fc_file_name, 'sub-(\d+)-FC', 'tokens');
    if isempty(tokens)
        continue;
    end
    fc_id = tokens{1}{1}; 
    fc_id_padded = sprintf('%05d', str2double(fc_id)); % e.g., '00009'
    subID = fc_id_padded;
    FCfiledir = fullfile(FCdirpath, fc_file_name);
    flagfc = 1;
    
    % Match corresponding SC file
    for ij = 1:length(SCfiles)
        if SCfiles(ij).isdir
            continue;
    end
        sc_file_name = SCfiles(ij).name;
        sc_tokens = regexp(sc_file_name, '^(\d+)_StructuralConnectome', 'tokens');
        if isempty(sc_tokens)
            continue;
        end
        sc_id = sc_tokens{1}{1}; 
        if strcmp(fc_id_padded, sc_id)
            SCfiledir = fullfile(SCdirpath, sc_file_name);
            flagsc = 1;
            break
        end
    end
    
    % Check if both modalities are successfully matched
    if flagfc * flagsc * flagmn == 1
        fc = importdata(FCfiledir);
        W = importdata(SCfiledir);
        Sc = W; Sc(Sc>0) = 1;              % Structural binary matrix
        N = size(W,1);                     % Node number (e.g., 400 for Schaefer)
        
        disp(['Processing Subject ID: ', subID]);
        disp('Multimodal data matched successfully.');
        numflag = numflag + 1;
        fc(fc<0) = 0;                      % Exclude negative functional connections
    else
        continue
    end
    
    % Load nodal coordinates
    load coor.mat 
    % coor represents x,y,z node coordinates | N x 3 matrix
  
    % ---------------------------------------------------------------------
    % 3. Graph Harmonic Model & Laplacian Eigenmodes
    % ---------------------------------------------------------------------
    A = W; A = -A;
    for i = 1:N
        A(i,i) = -sum(A(i,:));
    end
    B = A / max(eig(A));                   % Structural Laplacian matrix
    [BEC, BE] = eig(B);                    % BEC-eigenvectors, BE-eigenvalues

    [a, b] = sort(diag(BE));
    BE_sort = diag(a);
    BEC_sort = BEC(:, b);

    % Diffusion map algorithm on functional network
    mappedX = diffusion_maps(fc, 10, 0.5, 1);   
    x = mappedX(:, 1);                     % Node weights for first eigenvector
    xx = x * -1;                           % Reverse weights: positive = top of hierarchy
    allxx(:, numflag) = xx;
    
    % ---------------------------------------------------------------------
    % 4. Prediction Model (S-F Coupling Calculation)
    % ---------------------------------------------------------------------
    % A multilinear model is used to predict the functional connection profile
    % of every node based on the 2nd-14th eigenmodes of structural Laplacian network.
    rsq = zeros(N, 1); % Node-wise vector
    
    % Define FC response (y) and SC predictors (x)
    for jj = 1:N    
        y = fc(:, jj);   
        x1 = zscore(BEC_sort(:, 2:14));    % Standardize predictors of 2nd-14th eigenmodes    
        
        % Fit multiple regression (OLS, main effects only), excluding self-connections 
        lm = fitlm(x1, y, 'Exclude', jj);
        rsq(jj) = lm.Rsquared.Ordinary;    % Record R-squared for node jj
    end
    rsq1 = sqrt(real(rsq));                % Calculate S-F coupling value
    allrsq(:, numflag) = rsq1;

end % End of subject loop

% -------------------------------------------------------------------------
% 5. Visualization & Statistics
% -------------------------------------------------------------------------
figure(3)
rsq1 = BEC_sort(:,1);
rsq1 = rsq1';
x = coor(:,1); y = coor(:,2); z = coor(:,3);
N1 = ceil(max(rsq1)*10000) - floor(min(rsq1)*10000);
S = 200 * ones(size(x));

% Color mapping (Red-Blue)
mycolorpoint = [[255 135 158];[255 191 204];[250 235 214];[239 255 255];[64 104 224];[20 48 135]]/255;
mycolorpoint = flipud(mycolorpoint); 
mycolorposition = [1 1000 round(N1/2)-100 round(N1/2)+100 N1-1000 N1];
mycolormap_r = interp1(mycolorposition, mycolorpoint(:,1), 1:N1, 'linear', 'extrap');
mycolormap_g = interp1(mycolorposition, mycolorpoint(:,2), 1:N1, 'linear', 'extrap');
mycolormap_b = interp1(mycolorposition, mycolorpoint(:,3), 1:N1, 'linear', 'extrap');
C = [mycolormap_r', mycolormap_g', mycolormap_b']; 

colormap(C);
scatter3(x, y, z, S, C(round((rsq1-min(rsq1))*10000)+1,:), 'filled')
axis off; colorbar('Ticks',[])
figure(6); imagesc(rsq1);

% -------------------------------------------------------------------------
% Figure: R values anticorrelated with functional gradient
% -------------------------------------------------------------------------
figure;
allrsq(:, numflag+1:end) = [];
allxx(:, numflag+1:end) = [];
y1 = mean(allrsq, 2);
rsq1 = y1;
xx = mean(allxx, 2);

[rho, pval] = corr(xx, y1); 
lm = fitlm(xx, y1);
xhat = linspace(min(xx), max(xx), 100);
yhat = lm.Coefficients.Estimate(1) + (lm.Coefficients.Estimate(2) * xhat);

plot(xx, y1, '.', 'Markersize', 10); hold on
plot(xhat, yhat, 'LineWidth', 2)
title(['rho = ' num2str(rho) ', p = ' num2str(pval)]);
axis square
xlabel('Gradient'); ylabel('S-F Coupling (R)');

% -------------------------------------------------------------------------
% Boxplot of R across functional networks
% -------------------------------------------------------------------------
load lab.mat
load('rsn_mapping.mat');

R_name = rsn_names;                       % Total region names 
R_n = length(R_name);                     % Total region number
Rn_Ind = lab;                             % Network index for each node
Rn_name = cell(length(rsq1), 1);          % Network name for each node
R_R = zeros(R_n, 1);                      % Median R
R_A = zeros(R_n, 1);                      % Average R

for i = 1:R_n
    R_R(i) = median(rsq1(Rn_Ind==i));
    R_A(i) = mean(rsq1(Rn_Ind==i));
    Rn_name(Rn_Ind==i) = R_name(i);
end

[temp, I] = sort(R_R, 'descend');         % Sort networks by R median
a = rsq1; b = Rn_name;
c = R_name(I);                            % Sorted network names                    
Color = [219 2 10; 231 95 27; 238 146 43; 246 191 65; 246 236 84; 202 222 169; 147 205 137; 76 177 99]/255;
Color = flipud(Color);

figure;
valid_idx = ~cellfun(@isempty, b);
boxplot(a(valid_idx), b(valid_idx), 'Orientation', 'horizontal', 'GroupOrder', c, ...
    'Colors', Color, 'BoxStyle', 'filled', 'Whisker', 0, 'OutlierSize', 2.5, ...
    'Symbol', '.', 'MedianStyle', 'target');
title('Structural-Functional Coupling across Networks');
   
% -------------------------------------------------------------------------
% 6. Save output variables to CSV
% -------------------------------------------------------------------------
csv_file_path = 'Output_SFC_Values.csv';  % Specify output CSV file path
dlmwrite(csv_file_path, allrsq);
disp('---------------------------------------------------');
disp('Analysis completed. S-F coupling data successfully saved to CSV.');
disp('---------------------------------------------------');