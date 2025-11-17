% Script: 04_run_OptW.m
%
% Purpose: Optimizes the weight vector W for a composite risk score 
%          to maximize the Cox proportional hazards regression coefficient (beta).
%          A non-linear constraint is applied to ensure the Hazard Ratio (HR) 
%          confidence interval (CI) width remains below a threshold (100).
%
% Dependencies:
%   - Cox_Opt.m (Objective function)
%   - HRConfidenceConstraint.m (Non-linear constraint function)
%   - ../data/input/TCGA.csv (Input data file)
%
% Outputs:
%   - W_opt: Optimized weight vector
%   - Comparison of HRs (Optimized, Mean, Individual)
%
% ---

clear all
clc

% --- Configuration ---
% Ensure the current directory is the 'main' directory or the functions are in the path.
% We assume Cox_Opt.m and HRConfidenceConstraint.m are in the same folder.
INPUT_DATA_FILE = '../data/input/TCGA.csv';

% Optimization Parameters
CI_DIFFERENCE_THRESHOLD = 100; % Constraint: HR CI upper - HR CI lower <= 100
MAX_ITERATIONS = 1000;
MAX_FUNC_EVALS = 5000;
Z_SCORE_95CI = 1.96; % For 95% CI calculation

% --- Load Data ---
Mydata = importdata(INPUT_DATA_FILE);

% Survival Data
T = Mydata.data(:, 1);        % Survival/Follow-up time
status = Mydata.data(:, 2);   % 0=Event (Death), 1=Censor (Alive)
censor = status;              % 'Censoring' parameter expects 1 for censored observation

% Target Variables (Genes/Scores)
% ASSUMPTION: The target data (e.g., 6 genes/scores) starts at column 12490
% and spans 6 columns.
START_COL = 12490;
NUM_GENES = 6;
TargetData = Mydata.data(:, START_COL : (START_COL + NUM_GENES - 1));

% --- Optimization Setup ---

% Objective function: Maximize beta (i.e., Minimize -beta)
objective = @(W) -Cox_Opt(TargetData, W, T, censor);

% Initial point (W0): Equal weights sum to 1
W0 = ones(NUM_GENES, 1) / NUM_GENES; 

% Linear Equality Constraint: sum(W) = 1
Aeq = ones(1, NUM_GENES);
beq = 1;

% Bounds: W >= 0 and W <= 1
lb = zeros(NUM_GENES, 1); % Lower bound: weights must be non-negative
ub = ones(NUM_GENES, 1);  % Upper bound: weights must be <= 1

% Non-linear constraint function handle: HR Confidence Interval Constraint
nonlcon = @(W) HRConfidenceConstraint(W, TargetData, T, censor, CI_DIFFERENCE_THRESHOLD, Z_SCORE_95CI);

% Solver Options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxFunctionEvaluations', MAX_FUNC_EVALS, 'MaxIterations', MAX_ITERATIONS);

% Call fmincon
[W_opt, fval_opt] = fmincon(objective, W0, [], [], Aeq, beq, lb, ub, nonlcon, options);

% --- Validation and Results ---

% 1. Optimized Model Results
NewData = TargetData * W_opt;
[b_opt, ~, ~, stats_opt] = coxphfit(NewData, T, 'Censoring', censor); 
HR_opt = exp(b_opt);

% Calculate final HR CI for validation
if ~isempty(b_opt) && ~any(isnan(b_opt)) && stats_opt.se > 0
    b_lower_opt = b_opt - Z_SCORE_95CI * stats_opt.se;
    b_upper_opt = b_opt + Z_SCORE_95CI * stats_opt.se;
    hr_ci_opt_lower = exp(b_lower_opt);
    hr_ci_opt_upper = exp(b_upper_opt);
    hr_ci_diff_opt = hr_ci_opt_upper - hr_ci_opt_lower;
    
    disp(' ');
    disp('--- Optimization Results ---');
    disp(['Optimal HR (exp(beta)): ', num2str(HR_opt)]);
    disp(['Optimal 95% CI: [', num2str(hr_ci_opt_lower), ', ', num2str(hr_ci_opt_upper), ']']);
    disp(['Optimal CI Width: ', num2str(hr_ci_diff_opt)]);
    
    % Save Optimal Weights
    T_weights = array2table(W_opt, 'VariableNames', {'Weight'});
    writetable(T_weights, '../results/tables/OptW_weights.csv');
    disp('Optimal weights saved to ../results/tables/OptW_weights.csv');

else
    disp('--- Optimization Results ---');
    disp('Failed to calculate valid HR confidence interval for the optimal weights.');
end

% 2. Comparison: Average Score Model
MeanData = mean(TargetData, 2);
[b_mean] = coxphfit(MeanData, T, 'Censoring', censor);
HR_Mean = exp(b_mean);

% 3. Comparison: Individual Gene Models
HR_Each = zeros(1, NUM_GENES);
for i = 1:NUM_GENES
    [b] = coxphfit(TargetData(:, i), T, 'Censoring', censor);
    HR_Each(i) = exp(b);
end

% --- Display Comparisons ---
disp(' ');
disp('--- Comparison ---');
disp(['Optimal HR: ', num2str(HR_opt)]);
disp(['Average Score HR: ', num2str(HR_Mean)]);
disp('Individual Gene HRs:');
disp(HR_Each);

% Note on fval: fval_opt contains -b_opt, so exp(-fval_opt) should equal HR_opt.