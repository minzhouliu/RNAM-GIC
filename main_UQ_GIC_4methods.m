% Comparision of the calculation times of the four method for uncertainty 
%    quantification (UQ) of geomagnetically induced currents (GIC) in power grids:
%       1. Lehtinen-Pirjola method (LP method)
%       2. Nodal Admittance Matrix method (NAM method)
%       3. Bus Admittance Matrix method (BAM method)
%       4. Reduced Nodal Admittance Matrix method (RNAM method)
%
%
% Supplementary material for the letter submission:
%     Min-zhou Liu, Yan-zhao Xie, Yi-fan Yang, Riccardo Trinchero, and Igor S. Stievano, "Reduced 
%     Nodal Admittance Matrix Method for Probabilistic GIC Analysis in Power Grids," 2022.

clear;
clc;
close all;

%% Simulation configuration

% select a power grid test case
case_name = 'EPRI-21';
% case_name = 'IEEE_118-GMD';

% # of Monte Carlo samples for the uncertainty quantification of GICs
% n_MC = 1;
% n_MC = 10;
% n_MC = 100;
n_MC = 1e3;
% n_MC = 2e3;

% select the factorization method of the design matrix
factor_function = 1;  % sparse matrix factorization using MATLAB built-in function "decomposition"
% factor_function = 2;  % sparse matrix factorization using built-in function "lu"/"chol"
% factor_function = 3;  % full matrix factorization using built-in function "decomposition"
% factor_function = 4;  % full matrix factorization using built-in function "lu"/"chol"

% matrix factorization method
factor_method_LP = 'lu';       % LU factorization for LP method
factor_method_NAM = 'chol';    % Cholesky factorization for NAM method
factor_method_BAM = 'lu';      % LU factorization for BAM method
factor_method_RNAM = 'chol';   % Cholesky factorization for NAM method

% # of time instants
% n_t = 1;
% n_t = 100;
% n_t = 1440;
% n_t = 1440*12;
n_t = 1440*60;   % simulate GIC of 1-day, using 1-sec E-field data

% Seed used for random number generation to reproduce the results
rng('default');
% rng(2022);
% rng(2050);

% Large number for the virtual grounded branch's resistance of bus nodes in the LP method
M_big = 1e10;


%% Load the power grid matrices for GIC calculation

load(['.\data\Grid_matrices_', case_name, '.mat'],...
    'Y', 'n_g', 'n_b', 'Mat_incident', 'sub_no_remove', 'Mat_adjacency_bs', 'Rg');
% comments for the power grid matrices:
% Y is the full-node admittance matrix (including the admittance of grounded branches) in the NAM method
%       The arrangement of its blocks is Y = [Y_gg, Ygb; Y_bg, Y_bb];
% n_g is the # of substations in the DC circuit (The substations without grounded branch are not included, e.g.
%       substation 1 and 7 in the Benchmark EPRI-21 test case)
% n_b is the # of buses in the DC circuit (The generator bus or load bus connected 
%        to the delta winding side of the transformer is not included)
% Mat_incident is the incident matrix used to calculate the current injection at the bus nodes
%
% The following matrices or vectors are mainly used in the BAM method 
% sub_no_remove is a vector including the index number of the substations removed from the DC circuit
% Rg is a vector of the substation grounding resistance  (all substations are included)

Y_gg = Y(1:n_g, 1:n_g);     % block Y_gg
Y_bg = Y(n_g+1:end, 1:n_g);      % block Y_bg
Y_bb = Y(n_g+1:end, n_g+1:end);     % block Y_bb

n_g0 = length(Rg);      % # of all substations (n_g0 >= n_g)
sub_no_remain = setdiff(1:n_g0, sub_no_remove)';
Y_g = sparse(1:n_g, 1:n_g, 1./Rg(sub_no_remain));


%% Make the nodal current injections

% Load the GeoElectric Field (GEF) time series during the GMD event on August 15, 2015
% Column 1 is the northward GEF, and column 2 is the eastward GEF
load('.\data\GEF_YKC20150815_1sec.mat', 'GEF_ts');

% calculate current injections at bus nodes J_b using GEF time series 
J_b = Mat_incident * GEF_ts(:, 1:n_t);
% J_b = Mat_incident * GEF_ts(:, 1:60:end);
% synthetic current injections at bus nodes, and the random component is used 
%   to characterize spatial non-uniformity of GEF in a simple manner
% J_b = Mat_incident * GEF_ts(:, 1:n_t);   J_b = max(J_b, [], 'all') * 0.1 * randn(n_b, n_t);
% J_b = Mat_incident * GEF_ts(:, 1:n_t) + randn(n_b, n_t);
% J_b = (Mat_incident * [0;1]) .* ones(n_b, n_t);
% J_b = (Mat_incident * [1;0]) .* randn(n_b, n_t);
% J_b = (Mat_incident * [0;1]) .* randn(n_b, n_t);
% J_b = randn(n_b, n_t);


% Full-node current injections
J_full_node = [zeros(n_g, n_t); J_b];

%% Calulate GIC using four methods: LP, NAM, BAM & RNAM

% generate Monte Carlo samples of substation grounding resistances
Rg_samples = exp(randn(length(Rg), n_MC) * log(2)*0.5 + log(Rg));
% Verify GIC results using base substation grounding resistances in the references
% Rg_samples = Rg .* ones(n_g0, n_MC);

% simulation time of LP (column 1), NAM (column 2), BAM (column 3) & RNAM (column 4)
t_GIC = zeros(n_MC, 4);

%% GIC Calculation Method 1: LP

% network admittance matrix of LP method (Y_LP doesn't include grounded branches)
Y_LP = Y - diag([1./Rg(sub_no_remain); sparse(n_b,1)]);

for k = 1:n_MC
    
    % Rg_k is the k-th sample of substation grounding resistance vector
    Rg_k = Rg_samples(:,k);
    % update the grounding impedance matrix Z_full_node using new Rg_k
    Z_full_node = sparse(1:n_g+n_b, 1:n_g+n_b, [Rg_k(sub_no_remain); ones(n_b, 1)*M_big]);
    % update the design matrix of LP method
    Design_Mat_LP = speye(n_g+n_b) + Y_LP * Z_full_node;
    
    % solve the substation grounding GICs using the specified matrix factorization method
    tic
    if factor_function == 1
        Design_Mat_LP_factor = decomposition(Design_Mat_LP, factor_method_LP);
        GIC_LP = Design_Mat_LP_factor \ J_full_node;
    elseif factor_function == 2
        [L1, U1, P1, Q1, R1] = lu (Design_Mat_LP);
        GIC_LP = Q1 * (U1 \ (L1 \ (P1 * (R1 \ J_full_node))));
    elseif factor_function == 3
        Design_Mat_LP_factor = decomposition(full(Design_Mat_LP), factor_method_LP);
        GIC_LP = Design_Mat_LP_factor \ J_full_node;
    elseif factor_function == 4
        [L1, U1] = lu(full(Design_Mat_LP));
        GIC_LP = U1\(L1\J_full_node);
    end
    t_GIC(k,1) = toc;
    
end   
GIC_LP = GIC_LP(1:n_g, :);

%% GIC Calculation Method 2: NAM

for k = 1:n_MC
    
    % Rg_k is the k-th sample of substation grounding resistance vector
    Rg_k = Rg_samples(:,k);
    Design_Mat_NAM = Y;
    % update Y_gg using new Rg_k
    Design_Mat_NAM(1:n_g, 1:n_g) = Y_gg + sparse(1:n_g, 1:n_g, 1./Rg_k(sub_no_remain) - 1./Rg(sub_no_remain));
    % update the transformation matrix of NAM method
    Transform_Mat_NAM = sparse(1:n_g, 1:n_g, 1./Rg_k(sub_no_remain), n_g, n_g+n_b);
    
    % solve the substation grounding GICs using the specified matrix factorization method
    tic
    if factor_function == 1    
    %     Design_Mat_NAM_factor = decomposition(Design_Mat_NAM, factor_method_NAM);
        Design_Mat_NAM_factor = decomposition(Design_Mat_NAM, factor_method_NAM, 'upper');
        GIC_NAM = Transform_Mat_NAM  * (Design_Mat_NAM_factor\J_full_node);
    elseif factor_function == 2
        [L2, g2, P2] = chol (Design_Mat_NAM, 'lower') ;
        GIC_NAM = Transform_Mat_NAM * (P2 * (L2' \ (L2 \ (P2' * J_full_node))));
    elseif factor_function == 3
        Design_Mat_NAM_factor = decomposition(full(Design_Mat_NAM), factor_method_NAM, 'upper');
        GIC_NAM = Transform_Mat_NAM  * (Design_Mat_NAM_factor\J_full_node);
    elseif factor_function == 4
        R3 = chol(full(Design_Mat_NAM));
        GIC_NAM = Transform_Mat_NAM  * (R3\(R3'\J_full_node));
    end
    t_GIC(k,2) = toc;
    
end


%% GIC Calculation Method 3: BAM

% make matrices for BAM method
Ybn_BAM = diag(sum(-Y_bg, 2));
Ybb_BAM = Y_bb - diag(diag(Y_bb));
Ybb_BAM = Ybb_BAM - diag(sum(Ybb_BAM, 1));
Ybs_BAM = Ybb_BAM * Mat_adjacency_bs;
Ybs_BAM = Ybs_BAM(:, sum(Mat_adjacency_bs.*(1:n_g0),2));

for k = 1:n_MC
    
    % Rg_k is the k-th sample of substation grounding resistance vector
    Rg_k = Rg_samples(:,k);
    % update the matrix Ygb of BAM method using Rg_k
    Ygb_BAM = sparse(1:n_b, 1:n_b, Mat_adjacency_bs * (1./Rg_k));
    % update the design matrix and transformation matrix of BAM method
    Design_Mat_BAM = Ybn_BAM * Ygb_BAM + Ybs_BAM * Ybn_BAM + Ybb_BAM * Ygb_BAM;
    Transform_Mat_BAM = Mat_adjacency_bs(:,sub_no_remain)' * Ybn_BAM * Ygb_BAM;
    
    % solve the substation grounding GICs using the specified matrix factorization method
    tic
    if factor_function == 1
        Design_Mat_BAM_factor = decomposition(Design_Mat_BAM, factor_method_BAM);
        GIC_BAM = Transform_Mat_BAM * (Design_Mat_BAM_factor\J_b);
    elseif factor_function == 2
        [L3, U3, P3, Q3, R3] = lu (Design_Mat_BAM);
        GIC_BAM = Transform_Mat_BAM * ( Q3 * (U3 \ (L3 \ (P3 * (R3 \ J_b)))) );
    elseif factor_function == 3
        Design_Mat_BAM_factor = decomposition(full(Design_Mat_BAM), factor_method_BAM);
        GIC_BAM = Transform_Mat_BAM * (Design_Mat_BAM_factor\J_b);
    elseif factor_function == 4
        [L3, U3] = lu(full(Design_Mat_BAM));
        GIC_BAM = Transform_Mat_BAM * (U3\(L3\J_b));
    end
    t_GIC(k,3) = toc;
end


%% GIC Calculation Method 4: RNAM

for k = 1:n_MC
    
    % Rg_k is the k-th sample of substation grounding resistance vector
    Rg_k = Rg_samples(:,k);
    % the inverse of updated matrix block Y_gg
    Y_gg_inv = sparse(1:n_g, 1:n_g, 1./ (diag(Y_gg) + 1./Rg_k(sub_no_remain) - 1./Rg(sub_no_remain)) );
    % update the design matrix and transformation matrix of BAM method
    Design_Mat_RNAM = Y_bb - Y_bg * Y_gg_inv * Y_bg';
    Transform_Mat_RNAM = - sparse(1:n_g, 1:n_g, 1./Rg_k(sub_no_remain)) * Y_gg_inv * Y_bg';
    
    % solve the substation grounding GICs using the specified matrix factorization method
    tic
    if factor_function == 1
        Design_Mat_RNAM_factor = decomposition(Design_Mat_RNAM, factor_method_RNAM, 'upper');
        GIC_RNAM = Transform_Mat_RNAM * (Design_Mat_RNAM_factor \ J_b);
    elseif factor_function == 2
        [L4, g4, P4] = chol (Design_Mat_RNAM, 'lower') ;
        GIC_RNAM = Transform_Mat_RNAM * ( P4 * (L4' \ (L4 \ (P4' * J_b))) );
    elseif factor_function == 3
        Design_Mat_RNAM_factor = decomposition(full(Design_Mat_RNAM), factor_method_RNAM, 'upper');
        GIC_RNAM = Transform_Mat_RNAM * (Design_Mat_RNAM_factor \ J_b);
    elseif factor_function == 4
        R4 = chol(full(Design_Mat_RNAM), 'upper');
        GIC_RNAM = Transform_Mat_RNAM * (R4\(R4'\J_b));
    end
    t_GIC(k,4) = toc;

end

%% Post-Processing: Compare the GIC results of four methods
if n_g <= 10
    fprintf('Comparison of substation grounding GIC (A) at the 1st time instant :\n%15s %10s %10s %10s %10s\n',...
        'Substation #', 'LP', 'NAM', 'BAM', 'RNAM');
    fprintf('%15d %10.2f %10.2f %10.2f %10.2f\n',...
        full([sub_no_remain, GIC_NAM(:,1), GIC_LP(:,1), GIC_BAM(:,1), GIC_RNAM(:,1)]'));
else
    fprintf('Comparison of substation grounding GIC (A) at the 1st time instant (partical results):\n%15s %10s %10s %10s %10s\n',...
        'Substation #', 'LP', 'NAM', 'BAM', 'RNAM');
    fprintf('%15d %10.2f %10.2f %10.2f %10.2f\n',...
        full([sub_no_remain(1:10), GIC_NAM(1:10,1), GIC_LP(1:10,1), GIC_BAM(1:10,1), GIC_RNAM(1:10,1)]'));
end

figure;
plot(GIC_LP(:,1), 'LineWidth', 1.5);
hold on;
plot(GIC_NAM(:,1), '-.', 'LineWidth', 1.5);
plot(GIC_BAM(:,1), 'o', 'LineWidth', 1.5);
plot(GIC_RNAM(:,1), '--', 'LineWidth', 1.5);
xlabel('Substation #');
ylabel('Substation Grounding GIC (A)');
legend({'LP', 'NAM', 'BAM', 'RNAM'}, 'Location','best')
set(gca, 'XTick', 1:n_g, 'XTickLabel', sub_no_remain);
title('GIC at the 1st time instant')
xlim([0 min(25, n_g+1)])


%% Compare the calculation time of four methods
fprintf('Comparison of calculation times(s) by using n_MC = %d Monte Carlo samples :\n%10s %10s %10s %10s\n',...
    n_MC, 'LP', 'NAM', 'BAM', 'RNAM')
fprintf('%10.2f %10.2f %10.2f %10.2f\n',...
    sum(t_GIC(:,1)), sum(t_GIC(:,2)), sum(t_GIC(:,3)), sum(t_GIC(:,4)));

fprintf('Estimated results:\nComparison of calculation times(s) if # of Monte Carlo n_MC = 1000:\n%10s %10s %10s %10s\n',...
    'LP', 'NAM', 'BAM', 'RNAM')
fprintf('%10.2f %10.2f %10.2f %10.2f\n',...
    mean(t_GIC(:,1))*1e3, mean(t_GIC(:,2))*1e3, mean(t_GIC(:,3))*1e3, mean(t_GIC(:,4))*1e3);

% fprintf('Estimated results: Comparison of calculation times(s) if # of Monte Carlo n_MC = 2000:\n%10s %10s %10s %10s\n', 'LP', 'NAM', 'BAM', 'RNAM')
% fprintf('%10.2f %10.2f %10.2f %10.2f\n', mean(t_GIC(:,1))*2e3, mean(t_GIC(:,2))*2e3, mean(t_GIC(:,3))*2e3, mean(t_GIC(:,4))*2e3);

% Probability distribution of calculation time for four methods: histogram
figure;
histogram(t_GIC(:,1)*1e3)
hold on;
histogram(t_GIC(:,2)*1e3)
histogram(t_GIC(:,3)*1e3)
histogram(t_GIC(:,4)*1e3)
xlabel('Calculation time (ms)');
ylabel('Counts');
legend('LP', 'NAM', 'BAM', 'RNAM', 'Location','best')
set(gca, 'XScale', 'log')
grid on;

% Probability distribution of calculation time for four methods: boxplot
if n_MC > 1
    figure;
    boxplot(t_GIC*1e3);
    set(gca, 'XTickLabel', {'LP', 'NAM', 'BAM', 'RNAM'})
    xlabel('GIC calculation method')
    ylabel('Calculation time (ms)')
%     set(gca, 'YScale', 'log')
%     grid on;
end
