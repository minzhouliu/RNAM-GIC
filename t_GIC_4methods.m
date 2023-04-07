% Test and verify the results of four GIC calculation methods:
%       1. Lehtinen-Pirjola method (LP method)
%       2. Nodal Admittance Matrix method (NAM method)
%       3. Bus Admittance Matrix method (BAM method)
%       4. Reduced Nodal Admittance Matrix method (RNAM method)
%   and the results of the four methods are consistent with the public literature
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

% Large number for the virtual grounded branch's resistance of bus nodes in the LP method
M_big = 1e10;

% 1V/km uniform GeoElectric Field (GEF)
% E = [1;0];      % 1V/km northward GEF
E = [0;1];      % 1V/km eastward GEF

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

% calculate current injections at bus nodes J_b using GEF
J_b = Mat_incident * E;

% Full-node current injections
J_full_node = [zeros(n_g, size(J_b, 2)); J_b];

%% Load the results in the public reference to verify the code

% GIC results from the public literature:
% R. Horton, D. Boteler, T. J. Overbye, R. Pirjola, and R. C. Dugan, "A
% test case for the calculation of geomagnetically induced currents," IEEE
% Trans. Power Del., vol. 27, no. 4, pp. 2368-2373, 2012.

if strcmp(case_name, 'EPRI-21')
    if all(E==[1;0])
        GIC_ref = [0; 115.63; 139.85; 19.98; -279.08; -57.29; 0; 60.90];
    elseif all(E==[0;1])
        GIC_ref = [0; -189.29; -109.49; -124.58; -65.46; 354.52; 0; 134.30];
    end
end


%% GIC Calculation Method 1: LP

Z_full_node = sparse(1:n_g+n_b, 1:n_g+n_b, [Rg(sub_no_remain); ones(n_b, 1)*M_big]);
Y_LP = Y - diag([1./Rg(sub_no_remain); sparse(n_b,1)]);

Design_Mat_LP = speye(n_g+n_b) + Y_LP * Z_full_node;
Design_Mat_LP_factor = decomposition(Design_Mat_LP, 'lu');
GIC_LP = Design_Mat_LP_factor \ J_full_node;
GIC_LP = GIC_LP(1:n_g, :);

%% GIC Calculation Method 2: NAM

% make matrices for classical full-node NAM method
Design_Mat_NAM = Y;
Transform_Mat_NAM = [sparse(1:n_g, 1:n_g, 1./Rg(sub_no_remain)), sparse(n_g, n_b)];
Design_Mat_NAM_factor = decomposition(Design_Mat_NAM, 'chol', 'upper');
GIC_NAM = Transform_Mat_NAM  * (Design_Mat_NAM_factor\J_full_node);

%% GIC Calculation Method 3: BAM

% make matrices for BAM method
Ygb_BAM = sparse(1:n_b, 1:n_b, Mat_adjacency_bs * (1./Rg));
Ybn_BAM = diag(sum(-Y_bg, 2));
Ybb_BAM = Y_bb - diag(diag(Y_bb));
Ybb_BAM = Ybb_BAM - diag(sum(Ybb_BAM));
Ybs_BAM = Ybb_BAM * Mat_adjacency_bs;
Ybs_BAM = Ybs_BAM(:, sum(Mat_adjacency_bs.*(1:n_g0),2));

Design_Mat_BAM = Ybn_BAM * Ygb_BAM + Ybs_BAM*Ybn_BAM + Ybb_BAM*Ygb_BAM;

Transform_Mat_BAM = Mat_adjacency_bs(:,sub_no_remain)' * Ybn_BAM * Ygb_BAM;

Design_Mat_BAM_factor = decomposition(Design_Mat_BAM, 'lu');
GIC_BAM = Transform_Mat_BAM * (Design_Mat_BAM_factor\J_b);

% Compare the BAM Design matrix calculated with the results calculated
% through the code shared in the public literature:
% S. Marsal et al., "A new standalone tool for DC-equivalent network generation 
% and GIC calculation in power grids with multiple voltage levels," Space Weather, 
% vol. 20, no. 3, p. e2021SW002984, 2022.
% 
% if strcmp(case_name, 'EPRI-21')
%     load('.\data\Grid_matrices_EPRI_21_BAM_ref.mat', 'DesignMat_BAM_ref', 'bus_no_order');
%     fprintf('Norm of the difference between the BAM design matrix we calculated and the design matrix from the public literature: %.2e\n',...
%         norm(full(Design_Mat_BAM) - full(DesignMat_BAM_ref(bus_no_order, bus_no_order))));
%     
%     fig1 = figure;
%     subplot(121)
%     heatmap(Design_Mat_BAM);
%     title('BAM design matrix we calculated')
%     subplot(122)
%     heatmap(DesignMat_BAM_ref(bus_no_order, bus_no_order));
%     title('BAM design matrix from the public literature')
%     set(fig1, 'position', get(0, 'ScreenSize'))
% end


%% GIC Calculation Method 4: RNAM

Y_gg_inv = sparse(1:n_g, 1:n_g, 1./diag(Y_gg));

Design_Mat_RNAM = Y_bb - Y_bg * Y_gg_inv * Y_bg';
Transform_Mat_RNAM = - Y_g * Y_gg_inv * Y_bg';

Design_Mat_RNAM_factor = decomposition(Design_Mat_RNAM, 'chol', 'upper');
GIC_RNAM = Transform_Mat_RNAM * (Design_Mat_RNAM_factor \ J_b);


%% Compare the GIC results of four methods
if n_g <= 10
    fprintf('\nComparison of substation grounding GIC (A) at the 1st time instant :\n%15s %10s %10s %10s %10s\n',...
        'Substation #', 'LP', 'NAM', 'BAM', 'RNAM');
    fprintf('%15d %10.2f %10.2f %10.2f %10.2f\n',...
        full([(1:n_g)', GIC_NAM(:,1), GIC_LP(:,1), GIC_BAM(:,1), GIC_RNAM(:,1)]'));
else
    fprintf('\nComparison of substation grounding GIC (A) at the 1st time instant (partical results):\n%15s %10s %10s %10s %10s\n',...
        'Substation #', 'LP', 'NAM', 'BAM', 'RNAM');
    fprintf('%15d %10.2f %10.2f %10.2f %10.2f\n',...
        full([(1:10)', GIC_NAM(1:10,1), GIC_LP(1:10,1), GIC_BAM(1:10,1), GIC_RNAM(1:10,1)]'));
end

figure;
plot(GIC_LP(:,1), 'LineWidth', 1.5);
hold on;
plot(GIC_NAM(:,1), '-.', 'LineWidth', 1.5);
plot(GIC_BAM(:,1), 'o', 'LineWidth', 1.5);
plot(GIC_RNAM(:,1), '--', 'LineWidth', 1.5);
if strcmp(case_name, 'EPRI-21')
    scatter(1:n_g, GIC_ref(sub_no_remain), 100, '^', 'r');
end
xlabel('Substation #');
ylabel('Substation Grounding GIC (A)');
if strcmp(case_name, 'EPRI-21')
    legend({'LP', 'NAM', 'BAM', 'RNAM', 'Reference'}, 'Location','best')
else
    legend({'LP', 'NAM', 'BAM', 'RNAM'}, 'Location','best')
end
set(gca, 'XTick', 1:n_g, 'XTickLabel', sub_no_remain);
xlim([0 min(25, n_g+1)])


%% Compare the R condition number of four design matrices

% The calculation of condition number is time-consuming for large matrix (e.g. 2,000-node matrix)
if n_g + n_b < 2e3
%     fprintf('Comparison of reciprocal condition number:\n%10s %10s %10s %10s\n', 'LP', 'NAM', 'BAM', 'RNAM')
%     fprintf('%10.2e %10.2e %10.2e %10.2e\n',...
%         rcond(full(Design_Mat_LP)), rcond(full(Y)), rcond(full(Design_Mat_BAM)), rcond(full(Design_Mat_RNAM)));

    tic
    fprintf('Comparison of condition number:\n%20s %10s %10s %10s\n', 'LP', 'NAM', 'BAM', 'RNAM')
    fprintf('%20.2e %10.2e %10.2e %10.2e\n',...
        cond(full(Design_Mat_LP)), cond(full(Y)), cond(full(Design_Mat_BAM)), cond(full(Design_Mat_RNAM)));
    fprintf('%20.2f %10.2f %10.2f %10.2f\n',...
        cond(full(Design_Mat_LP)), cond(full(Y)), cond(full(Design_Mat_BAM)), cond(full(Design_Mat_RNAM)));
    toc
end

% fprintf('Comparison of condition number:\n%10s %10s %10s %10s\n', 'LP', 'NAM', 'BAM', 'RNAM')
% fprintf('%10.2e %10.2e %10.2e %10.2e\n',...
%     condest(Design_Mat_LP), condest(Y), condest(Design_Mat_BAM), condest(Design_Mat_RNAM));


%% Compare the storage of four design matrices

storage_LP = whos('Design_Mat_LP').bytes/1024;
storage_NAM = whos('Design_Mat_NAM').bytes/1024;
storage_BAM = whos('Design_Mat_BAM').bytes/1024;
storage_RNAM = whos('Design_Mat_RNAM').bytes/1024;
fprintf('Comparison of storage (kB) of the Sparse design matrices:\n%10s %10s %10s %10s\n', 'LP', 'NAM', 'BAM', 'RNAM')
% fprintf('%10.2e %10.2e %10.2e %10.2e\n',...
%         storage_LP, storage_NAM, storage_BAM, storage_RNAM);
fprintf('%10.2f %10.2f %10.2f %10.2f\n',...
        storage_LP, storage_NAM, storage_BAM, storage_RNAM);

full_Design_Mat_LP = full(Design_Mat_LP);
full_Design_Mat_NAM = full(Design_Mat_NAM);
full_Design_Mat_BAM = full(Design_Mat_BAM);
full_Design_Mat_RNAM = full(Design_Mat_RNAM);
storage_full_LP = whos('full_Design_Mat_LP').bytes/1024;
storage_full_NAM = whos('full_Design_Mat_NAM').bytes/1024;
storage_full_BAM = whos('full_Design_Mat_BAM').bytes/1024;
storage_full_RNAM = whos('full_Design_Mat_RNAM').bytes/1024;
fprintf('Comparison of storage (kB) of the Full design matrices:\n%10s %10s %10s %10s\n', 'LP', 'NAM', 'BAM', 'RNAM')
% fprintf('%10.2e %10.2e %10.2e %10.2e\n',...
%         storage_full_LP, storage_full_NAM, storage_full_BAM, storage_full_RNAM);
fprintf('%10.2f %10.2f %10.2f %10.2f\n',...
        storage_full_LP, storage_full_NAM, storage_full_BAM, storage_full_RNAM);
    
