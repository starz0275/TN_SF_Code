% =========================================================================
% Transcriptomic Association using Partial Least Squares (PLS) Regression
% with Spatial Permutation Testing (Spin-test)
% =========================================================================
clear; clc; close all;

% Add ENIGMA toolbox or relevant spin-test functions to path
% addpath(genpath('./ENIGMA_MacBrain_Toolbox')); 

% -------------------------------------------------------------------------
% 0. Environment Setup & Data Loading
% -------------------------------------------------------------------------
% Define relative paths (Please ensure these files are in the repository)
gene_path   = './Data_Sample/AHBA_Schaefer400_Gene_Expression.csv';
map_path    = './Data_Sample/Treatment_Effect_for_PLS_Tonly.csv';
coords_path = './Data_Sample/Schaefer400_Centroid_RAS.csv';
output_dir  = './Output_PLS_Results/';

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

fprintf('--- Loading Multimodal & Transcriptomic Data ---\n');

% Extract gene names from the header
fid = fopen(gene_path, 'r', 'n', 'UTF-8');
header_line = fgetl(fid); 
fclose(fid);
gene_names_raw = strsplit(header_line, ',');
gene_names = gene_names_raw(2:end); 

% Load numerical matrices (skipping header and first column labels)
X = dlmread(gene_path, ',', 1, 1);
Y_raw = dlmread(map_path);
Y = Y_raw(:, end); % Extract the neuroimaging effect size map (e.g., T-values)

% Z-score normalization for both transcriptomic and neuroimaging data
X_scaled = zscore(X);
Y_scaled = zscore(Y);
n_region = size(X_scaled, 1);
n_gene = size(X_scaled, 2);

% -------------------------------------------------------------------------
% 1. Multi-component PLS Analysis (Evaluating Variance Explained)
% -------------------------------------------------------------------------
n_comp = 10;
fprintf('--- Running %d-Component PLS Regression ---\n', n_comp);
[XL, YL, XS, YS, BETA, PCTVAR] = plsregress(X_scaled, Y_scaled, n_comp);

% Calculate observed correlation for PLS1
r_obs = corr(XS(:, 1), Y_scaled);
y_explained_each = PCTVAR(2, :) * 100;    % Variance explained by each component (%)
y_cumulative = cumsum(y_explained_each);  % Cumulative variance explained (%)

fprintf('Variance explained by PLS1: %.2f%%\n', y_explained_each(1));
fprintf('Cumulative variance explained by first 5 components: %.2f%%\n', y_cumulative(5));

% -------------------------------------------------------------------------
% 2. Full PLS Spin-test (Spatial Autocorrelation Correction)
% -------------------------------------------------------------------------
n_spin = 1000;
null_pls_all = zeros(n_spin, n_comp); 
p_spin_all = zeros(n_comp, 1);

fprintf('--- Executing Full PLS Spin-test (%d permutations) ---\n', n_spin);

% Step 2.1: Generate spatially constrained null models (Spin rotations)
% Load regional centroids to define spatial coordinates
coords_table = readtable(coords_path);
lh_centroid = [coords_table.R(1:200), coords_table.A(1:200), coords_table.S(1:200)];
rh_centroid = [coords_table.R(201:400), coords_table.A(201:400), coords_table.S(201:400)];

% Generate spin permutation indices (e.g., using ENIGMA's rotate_parcellation logic)
% perm_id matrix size: (400 regions x n_spin permutations)
perm_id = rotate_parcellation(lh_centroid, rh_centroid, n_spin);

% Step 2.2: Execute permutation loop
for i = 1:n_spin
    if mod(i, 100) == 0; fprintf('  Progress: %d / %d permutations\n', i, n_spin); end
    
    % Reorder the neuroimaging map using spatial spin indices
    Y_spin = Y_scaled(perm_id(:, i)); 
    
    % Re-run the PLS regression on the spun data
    [~, ~, XS_spin] = plsregress(X_scaled, Y_spin, n_comp);
    
    % Store the pseudo-correlation coefficients for each component
    for k = 1:n_comp
        null_pls_all(i, k) = corr(XS_spin(:, k), Y_spin);
    end
end

% Step 2.3: Calculate the spin-based P-values (P_spin)
for k = 1:n_comp
    r_obs_k = corr(XS(:, k), Y_scaled); % True observed correlation
    % Non-parametric P-value calculation
    p_spin_all(k) = (sum(abs(null_pls_all(:, k)) >= abs(r_obs_k)) + 1) / (n_spin + 1);
end

% Display spin-test summary
results_summary = table((1:n_comp)', y_explained_each', y_cumulative', p_spin_all, ...
    'VariableNames', {'Component', 'Individual_Var', 'Cumulative_Var', 'P_spin'});
disp(results_summary);

% -------------------------------------------------------------------------
% 3. Gene Stability Analysis (Bootstrap Resampling)
% -------------------------------------------------------------------------
n_boot = 1000;
boot_weights = zeros(n_gene, n_boot);
fprintf('--- Running Bootstrap Resampling for Gene Stability (n=%d) ---\n', n_boot);

% Extract the original PLS1 gene weights for sign alignment
gene_weight_real = XL(:, 1); 

parfor i = 1:n_boot % Parallel computing to accelerate
    idx = randsample(n_region, n_region, true); % Resampling with replacement
    [XLb, ~, ~, ~, ~, ~] = plsregress(X_scaled(idx, :), Y_scaled(idx), 1);
    
    current_weight = XLb(:, 1);
    
    % Procrustes rotation / Sign-flip correction: 
    % Ensures the bootstrap component directions align with the original PLS1
    if dot(current_weight, gene_weight_real) < 0
        current_weight = -current_weight;
    end
    
    boot_weights(:, i) = current_weight;
end

gene_weight = gene_weight_real;
gene_std = std(boot_weights, 0, 2);
gene_Z = gene_weight ./ gene_std; % Calculate Bootstrap Z-scores (Stability)

% -------------------------------------------------------------------------
% 4. Results Export
% -------------------------------------------------------------------------
res_table = table(gene_names', gene_weight, gene_Z, ...
                  'VariableNames', {'Gene', 'Weight', 'Z_score'});
res_sorted = sortrows(res_table, 'Z_score', 'descend');

% Export full results and top 5% driving genes
writetable(res_sorted, fullfile(output_dir, 'PLS1_Full_Gene_Results.csv'));
num_top = round(0.05 * n_gene);
writetable(res_sorted(1:num_top, :), fullfile(output_dir, 'PLS1_Top_5_Percent_Genes.csv'));

fprintf('--- Analysis Completed. Results saved to output directory. ---\n');

% -------------------------------------------------------------------------
% 5. Comprehensive Visualization
% -------------------------------------------------------------------------
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.9 0.4]);

% Subplot A: Variance Explained (10 Components)
subplot(1, 2, 1);
yyaxis left
bar(1:n_comp, y_explained_each, 'FaceColor', [0.5 0.2 0.5], 'FaceAlpha', 0.3);
ylabel('Individual Variance (%)');
yyaxis right
plot(1:n_comp, y_cumulative, '-o', 'Color', [0.5 0 0.5], 'LineWidth', 2);
ylabel('Cumulative Variance (%)');
title('A. Variance Explained across PLS Components');
xlabel('PLS Component'); 

% Subplot B: Spin-test Null Distribution (PLS1)
subplot(1, 2, 2);
r_obs_pls1 = corr(XS(:, 1), Y_scaled); 
histogram(null_pls_all(:, 1), 50, 'Normalization', 'pdf', 'FaceColor', [0.2 0.33 0.49], 'EdgeColor', 'w');
hold on;
line([r_obs_pls1 r_obs_pls1], get(gca, 'YLim'), 'Color', 'r', 'LineWidth', 3, 'LineStyle', '--');
title(sprintf('B. PLS1 Spin-test ($P_{spin}$  = %.4f)', p_spin_all(1)),'Interpreter', 'latex');
xlabel('Correlation Coefficient'); ylabel('Density');
legend('Null Distribution', sprintf('Observed r = %.3f', r_obs_pls1), 'Location', 'best');