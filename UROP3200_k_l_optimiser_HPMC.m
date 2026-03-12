% MAIN SCRIPT
disp('HPMC k_l optimisers');

% Optimised parameter from control (no polymer)
B_s = 2.0000e+03;
k_g = 1.0000e-04;
g = 1.7500;
zhat_mult = 0.2000;

% Experimental time range
exp_tt = 60 * [1, 2, 3,	4, 5, 10, 15, 20, 30, 45, 60, 90, 120, 150, 180];

% Range of polymer concentration that were tested
C_P = [0;1;2;5];

% Experimental bulk drug concentration 
exp_CC = [0.371, 0.456, 0.523, 0.554, 0.595, 0.722, 0.594, 0.498, 0.454, 0.438, 0.426, 0.406, 0.403, 0.400, 0.397;
            0.351, 0.452, 0.505, 0.543, 0.584, 0.694, 0.778, 0.832, 0.967, 0.828, 0.749, 0.659, 0.626, 0.600, 0.580;
            0.364, 0.465, 0.512, 0.557, 0.587, 0.701, 0.767, 0.817, 0.882, 0.926, 0.943, 0.946, 0.935, 0.954, 0.855;
            0.327, 0.433, 0.479, 0.520, 0.541, 0.674, 0.781, 0.810, 0.905, 0.939, 0.963, 0.965, 0.970, 0.985, 0.969];

% Other parameters to set
w_ccini = 303.37; % [=] mg
C_b_ini = exp_CC(:,1);
k_l_initial = [5, 5];

% HPMC_k_l_optimiser(w_ccini, C_b_ini, B_s, k_g, g, zhat_mult, C_P, exp_tt, exp_CC, k_l_initial)

% [model_CC, CC_full, tt_full] = UROP2100m1SCALED_final(w_ccini, C_b_ini, C_P, k_l_initial, B_s, k_g, g, zhat_mult, exp_tt, exp_CC, true);

% HPMC_k_l_optimiser(w_ccini, C_b_ini, B_s, k_g, g, zhat_mult, C_P, exp_tt, exp_CC, k_l_initial)

[k_l_opt, model_CC, CC_full, tt_full] = HPMC_k_l_optimiser(w_ccini, C_b_ini, B_s, k_g, g, zhat_mult, C_P, exp_tt, exp_CC, k_l_initial);

% Graphing
% Create figure with 2x2 subplots
figure('Position', [100, 100, 1200, 800]);

% Generate colors for each subplot
colors = lines(4);

% Loop through each C_P value (4 total)
for i = 1:length(C_P)
    % Create subplot (2 rows, 2 columns, position i)
    subplot(2, 2, i);
    hold on;
    
    % Plot 1: Full predicted curve (continuous line)
    plot(tt_full{i}/60, CC_full{i}, 'k-', 'LineWidth', 1.5, ...
         'DisplayName', 'Predicted curve');
    
    % Plot 2: Experimental data points (red x)
    plot(exp_tt/60, exp_CC(i, :), 'rx', 'LineWidth', 1, 'MarkerSize', 8, ...
         'DisplayName', 'Experimental');
    
    % Plot 3: Model data points at experimental times (blue x)
    plot(exp_tt/60, model_CC(i, :), 'bx', 'LineWidth', 1, 'MarkerSize', 8, ...
         'DisplayName', 'Model');
    
    % Customize the subplot
    xlabel('Time (min)');
    ylabel('Concentration (mg/ml)');
    
    % Create title based on C_P value
    if C_P(i) == 0
        title(sprintf('Control (C_P = 0 mg/L)'), 'FontWeight', 'bold');
    else
        title(sprintf('C_P = %.1f mg/L', C_P(i)), 'FontWeight', 'bold');
    end
    
    legend('show', 'Location', 'best');
    grid on;
    
    % Set axis limits if needed (adjust based on your data)
    % xlim([0, 180]);
    % ylim([0, max(CC_full{i}) * 1.1]);
    
    hold off;
end

% Add an overall title to the figure
sgtitle(sprintf('Cocrystal Dissolution - Polymer Concentration Effect\nOptimal k_l_{nuc} = %.4f, k_l_{growth} = %.4f', k_l_opt(1), k_l_opt(2)), ...
        'FontSize', 14, 'FontWeight', 'bold');

function [k_l_opt, model_CC, CC_full, tt_full] = HPMC_k_l_optimiser(w_ccini, C_b_ini, B_s, k_g, g, zhat_mult, C_P, exp_tt, exp_CC, k_l_initial)

polymer_idx = find(C_P > 0);

fprintf('\n=== 2-Parameter Optimization: k_l_nuc, k_l_growth ===\n');
fprintf('C_P=0 excluded from objective (polymer inhibition inactive for control)\n');
fprintf('Separate inhibition: nucleation (B_s) vs growth (G_s, G_b)\n\n');

% --- Step 1: 2D grid search ---
n_grid = 15;
k_l_nuc_vals = logspace(-3, 2, n_grid);
k_l_growth_vals = logspace(-3, 2, n_grid);
err_grid = inf(n_grid, n_grid);

fprintf('2D Grid search (%dx%d = %d evaluations)...\n', n_grid, n_grid, n_grid^2);
for i = 1:n_grid
    for j = 1:n_grid
        err_grid(i,j) = error_func([k_l_nuc_vals(i), k_l_growth_vals(j)]);
    end
end

[best_grid_err, min_idx] = min(err_grid(:));
[best_i, best_j] = ind2sub([n_grid, n_grid], min_idx);
best_k_l = [k_l_nuc_vals(best_i), k_l_growth_vals(best_j)];
best_err = best_grid_err;
fprintf('  Best grid: k_l_nuc=%.4f, k_l_growth=%.4f, SSE=%.6f\n', best_k_l(1), best_k_l(2), best_grid_err);

% Plot 2D SSE landscape
figure('Name', 'k_l 2D Optimization Landscape');
contourf(log10(k_l_growth_vals), log10(k_l_nuc_vals), log10(err_grid), 20);
colorbar;
xlabel('log_{10}(k_l_{growth})'); ylabel('log_{10}(k_l_{nuc})');
title('log_{10}(SSE) Landscape'); hold on;
plot(log10(best_k_l(2)), log10(best_k_l(1)), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
hold off;

% --- Step 2: Multi-start fminsearch ---
opts_fms = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxIter', 10000, ...
    'MaxFunEvals', 20000, 'Display', 'off');

% Collect top grid points as seeds
[~, sorted_idx] = sort(err_grid(:));
seeds = zeros(0, 2);
for s = 1:min(8, numel(sorted_idx))
    [si, sj] = ind2sub([n_grid, n_grid], sorted_idx(s));
    seeds(end+1, :) = [k_l_nuc_vals(si), k_l_growth_vals(sj)];
end

manual_seeds = [k_l_initial(:)'; 0.01 0.01; 0.1 0.1; 1 1; 5 5; 10 10; ...
                0.01 1; 1 0.01; 0.1 10; 10 0.1; 50 50; 0.5 5; 5 0.5; ...
                0.001 0.1; 0.1 0.001; 100 1; 1 100; 0.05 0.5; 0.5 0.05];
seeds = [seeds; manual_seeds];

fprintf('Running fminsearch from %d seeds...\n', size(seeds,1));
for s = 1:size(seeds,1)
    try
        [k_trial, e_trial] = fminsearch(@error_func_bounded, seeds(s,:), opts_fms);
        if e_trial < best_err && all(k_trial > 0)
            best_err = e_trial;
            best_k_l = k_trial;
            fprintf('  seed [%.3f, %.3f]: k_l_nuc=%.6f, k_l_growth=%.6f, SSE=%.6f *\n', ...
                seeds(s,1), seeds(s,2), k_trial(1), k_trial(2), e_trial);
        end
    catch
    end
end

k_l_opt = best_k_l;

% --- Final evaluation ---
[model_CC, CC_full, tt_full] = UROP2100m1SCALED_final(w_ccini, C_b_ini, C_P, k_l_opt, B_s, k_g, g, zhat_mult, exp_tt, exp_CC, true);

% R^2 for polymer cases
ss_res_p = sum(sum((exp_CC(polymer_idx,:) - model_CC(polymer_idx,:)).^2));
ss_tot_p = sum(sum((exp_CC(polymer_idx,:) - mean(exp_CC(polymer_idx,:), 'all')).^2));
R2_polymer = 1 - (ss_res_p / ss_tot_p);

% R^2 overall
ss_res_all = sum(sum((exp_CC - model_CC).^2));
ss_tot_all = sum(sum((exp_CC - mean(exp_CC, 'all')).^2));
R2_total = 1 - (ss_res_all / ss_tot_all);

fprintf('\n=== Results ===\n');
fprintf('Optimal k_l_nuc:    %.6f\n', k_l_opt(1));
fprintf('Optimal k_l_growth: %.6f\n', k_l_opt(2));
fprintf('R^2 (polymer cases, C_P>0): %.4f\n', R2_polymer);
fprintf('R^2 (all cases):            %.4f\n', R2_total);

fprintf('\nPer-curve R^2:\n');
for i = 1:length(C_P)
    ss_res_i = sum((exp_CC(i,:) - model_CC(i,:)).^2);
    ss_tot_i = sum((exp_CC(i,:) - mean(exp_CC(i,:))).^2);
    R2_i = 1 - (ss_res_i / ss_tot_i);
    fprintf('  C_P = %4.1f: R^2 = %.4f\n', C_P(i), R2_i);
end

fprintf('\nPolymer inhibition at optimal parameters:\n');
for i = 1:length(C_P)
    inh_nuc = C_P(i) / (k_l_opt(1) + C_P(i));
    inh_growth = C_P(i) / (k_l_opt(2) + C_P(i));
    fprintf('  C_P = %4.1f: nucleation %.1f%%, growth %.1f%%\n', C_P(i), inh_nuc*100, inh_growth*100);
end

    function err = error_func(k_l_val)
        if any(k_l_val <= 0)
            err = 1e10;
            return;
        end
        try
            [mc, ~, ~] = UROP2100m1SCALED_final(w_ccini, C_b_ini, C_P, k_l_val, B_s, k_g, g, zhat_mult, exp_tt, exp_CC, false);
            err = sum(sum((mc(polymer_idx,:) - exp_CC(polymer_idx,:)).^2));
        catch
            err = 1e10;
        end
    end

    function err = error_func_bounded(k_l_val)
        if any(k_l_val <= 0)
            err = 1e10 + sum(abs(k_l_val(k_l_val <= 0))) * 1e8;
            return;
        end
        err = error_func(k_l_val);
    end

end

function [pred_CC, CC_full, tt_full] = UROP2100m1SCALED_final(w_ccini, C_b_ini, C_P, k_l, B_s, k_g, g, zhat_mult, exp_tt, exp_CC, graphing)

% Check if inputs are vectors (4x1 column vectors)
    if numel(C_b_ini) > 1 || numel(C_P) > 1
        % Ensure both are column vectors
        C_b_ini = C_b_ini(:);  % Reshape to column vector
        C_P = C_P(:);           % Reshape to column vector
        
        % Check if they have the same length
        if length(C_b_ini) ~= length(C_P)
            error('C_b_ini and C_P must have the same length');
        end
        
        n_pairs = length(C_b_ini);  % This will be 4 in your case
        
        % Pre-allocate output matrix: n_pairs ¡Á 15
        % 15 is the number of experimental time points
        n_time_points = 15;
        pred_CC = zeros(n_pairs, n_time_points);
        CC_full = cell(n_pairs, 1);  % Cell array for full concentration curves
        tt_full = cell(n_pairs, 1);  % Also store time vectors

        for i = 1:n_pairs
            fprintf('Running simulation %d of %d: C_b_ini = %.3f, C_P = %.2f\n', ...
                    i, n_pairs, C_b_ini(i), C_P(i));
            
            % Call the function recursively with scalar inputs
            % Set graphing to false for individual runs to avoid multiple figures
            [pred_row, CC_row, tt_row] = UROP2100m1SCALED_final(w_ccini, C_b_ini(i), C_P(i), ...
                                               k_l, B_s, k_g, g, zhat_mult, exp_tt, exp_CC, graphing);
            
            % Store the result in the output matrix
            % pred_row should be a 1¡Á15 vector
            if size(pred_row, 1) > size(pred_row, 2)
                pred_row = pred_row';  % Ensure it's a row vector
            end
            pred_CC(i, :) = pred_row;

            % Store full curve (as column vectors)
            CC_full{i} = CC_row(:);  % Ensure column vector
            tt_full{i} = tt_row(:);  % Ensure column vector

        end
        % if graphing
        %     hold on;
        %     % Plot experimental data for reference
        %     plot(exp_tt/60, exp_CC, 'rx', 'LineWidth', 0.5, 'DisplayName', 'Experimental');
        % 
        %     xlabel('Time (min)');
        %     ylabel('Concentration (mg/ml)');
        %     title('Cocrystal Dissolution - Multiple C_b_{ini} and C_P Pairs');
        %     legend('show', 'Location', 'best');
        %     grid on;
        %     hold off;
        % end
        
        return;  % Exit the function after handling vector inputs
    end

% w_ccini = 303.37;    %initial mass of cocrystal
% C_b_ini = 0.371;
% C_P = 0;
% k_l = 0.1;           % polymer adsorption
% B_s = 6.16*10^2;
% k_g = 3.86*10^(-5);
% g = 4.02;   % try fitting this (usually between 1~2)
% graphing = true;

% we want to scale everything 
% hence we are going to make everything dimensionless

% Scale factors

t_char = 1; 

mu_s0_char = 10^7;
mu_s1_char = 10^9;
mu_s2_char = 10^10;
mu_s3_char = 10^9;
w_cc_char = 10^2;
mu_b0_char = 10^3;
mu_b1_char = 10^2;
mu_b2_char = 10^3;
mu_b3_char = 1;
C_db_char = 1;

% rho = crystal mass density
% k_v = crystal volume shape factor
% drug mass fraction of the stable drug form

% C_ds is constant

% parameters

V_t = 200; % (ml=cm^3) % V_t = V_s + V_b
r_ccini = 0.0075; % (cm)
rho = 1100; % (mg cm^-3)
f_dsf = 0.87;
f_dcc = 0.80;
k_v = 3.20*10^-4;

wsingle_ccini = rho * (4/3) * pi *r_ccini^3;

no_cc = w_ccini / wsingle_ccini;

K_sp = 6.60*10^(-18)*236269^2 * 118090; % (mM)^3 multiply with MW_API <- CBZ (236269mg/mol)and MW_CF <- succinic acid (118090mg/mol)

%B_s = 6.16*10^2 * 10;%6.16*10^2; % (mL^-1)
k_nub = 0; % (mL^-1 = cm^-3)
n_ub = 0;
%k_g = 3.86*10^(-5); % (cm s^-1)
%g = 4.02;   % try fitting this (usually between 1~2)
zhat = zhat_mult * 3.80*10^(-5); % (cm^4 mg^-1 s^-1)
AAD = 9.1; % (%)

Csol_d = 0.4; % (mg/ml) (assumed value with 5% ethanol); drug solubility

% Calculation for Concentration of surface region (C_s)
% Cstosat_d ^ p * Cstosat_cf ^ q = K_sp
% Cstosat_d/p = Cstosat_cf/q
% p = 2, q = 1
% Hence Cstosat_d = 2*Cstosat_cf
% (2*Cstosat_cf) ^ 2 * Cstosat_cf  = K_sp
% 4 * Cstosat_cf ^ 3 = K_sp

Cstosat_cf = (K_sp / 4)^(1/3);
Cstosat_d = 2*Cstosat_cf;

C_ds = 1.20; % concentration of drug on surface

if isscalar(k_l)
    k_l = [k_l, k_l];
end

params = struct(...
        't_char', t_char,  ...
        'V_t', V_t, 'r_ccini', r_ccini, 'rho', rho, ...
        'f_dsf', f_dsf, 'f_dcc', f_dcc, 'k_v', k_v, ...
        'wsingle_ccini', wsingle_ccini, 'no_cc', no_cc, ...
        'K_sp', K_sp, 'B_s', B_s, 'k_nub', k_nub, 'n_ub', n_ub, ...
        'k_g', k_g, 'g', g, 'zhat', zhat, 'AAD', AAD, ...
        'Cstosat_cf', Cstosat_cf, 'Cstosat_d', Cstosat_d, 'C_ds', C_ds, ...
        'Csol_d', Csol_d, 'k_l', k_l);

tf = 180 * 60;

t_scaled_final = tf / t_char;

%opt = odeset('AbsTol', 1e-3);
%opt = odeset('MinStep', t_scaled_final/100000);
%[tt_scaled, XX_scaled] = ode45(@(t, X) f(t, X, params), [0, t_scaled_final], [0;0;0;0;w_ccini/w_char;0;0;0;0;C_b_ini/C_char], opt);


try
    % Try ode15s first
    % options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [tt_scaled, XX_scaled] = ode15s(@(t, y) f(t, y, params), [60, t_scaled_final], [0;0;0;0;w_ccini/w_cc_char;0;0;0;0;C_b_ini/C_db_char]);
catch
    % If ode15s fails, try ode45
    warning('ode45 failed, switching to ode15s');
    % options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [tt_scaled, XX_scaled] = ode45(@(t, y) f(t, y, params), [60, t_scaled_final], [0;0;0;0;w_ccini/w_cc_char;0;0;0;0;C_b_ini/C_db_char]);
end


% [tt_scaled, XX_scaled] = ode15s(@(t, X) f(t, X, params), [0, t_scaled_final], [0;0;0;0;w_ccini/w_cc_char;0;0;0;0;C_b_ini/C_db_char]);

CC = XX_scaled(:,10) * C_db_char;

tt = tt_scaled * t_char;

hold on;

% if (graphing)
%     hold on;
%     plot(tt,CC, 'r-','DisplayName', 'Mean plot');
%     xlabel('time');
%     ylabel('Concentration in bulk (mg/ml)');
%     legend('show');
%     title('Concentration against time - k_l:', k_l);
% 
% 
%     plot(exp_tt,exp_CC,'bx','DisplayName', 'Actual datapoint');
%     hold on;
% 
% end

pred_CC = interp1(tt, CC, exp_tt, 'linear');

CC_full = XX_scaled(:,10) * C_db_char;

tt_full = tt_scaled * t_char;

% hold off

function dXdt_scaled = f(t,X,~)

        pp = odeset;  % Get the options structure

        % Unpack parameters
        t_char = params. t_char;
        V_t = params.V_t;
        rho = params.rho;
        f_dsf = params.f_dsf;
        f_dcc = params.f_dcc;
        k_v = params.k_v;
        wsingle_ccini = params.wsingle_ccini;
        no_cc = params.no_cc;
        B_s = params.B_s;
        k_nub = params.k_nub;
        n_ub = params.n_ub;
        k_g = params.k_g;
        g = params.g;
        zhat = params.zhat;
        C_ds = params.C_ds;
        Csol_d = params.Csol_d;
        k_l_nuc = params.k_l(1);
        k_l_growth = params.k_l(2);

        % moment equation for surface
        mu_s0_scaled = X(1,1);
        mu_s1_scaled = X(2,1);
        mu_s2_scaled = X(3,1);
        mu_s3_scaled = X(4,1);

        mu_s0 = mu_s0_scaled * mu_s0_char;
        mu_s1 = mu_s1_scaled * mu_s1_char;
        mu_s2 = mu_s2_scaled * mu_s2_char;
        mu_s3 = mu_s3_scaled * mu_s3_char;

        w_cc_scaled = X(5,1);

        w_cc = w_cc_scaled * w_cc_char;

        % moment equation for bulk
        mu_b0_scaled = X(6,1);
        mu_b1_scaled = X(7,1);
        mu_b2_scaled = X(8,1);
        mu_b3_scaled = X(9,1);

        mu_b0 = mu_b0_scaled * mu_b0_char;
        mu_b1 = mu_b1_scaled * mu_b1_char;
        mu_b2 = mu_b2_scaled * mu_b2_char;
        mu_b3 = mu_b3_scaled * mu_b3_char;

        C_db_scaled = X(10,1);

        C_db = C_db_scaled * C_db_char;

        S_b = C_db/Csol_d;

        S_s = C_ds/Csol_d;

        B_b = k_nub*(S_b-1)^n_ub;

        polymer_inh_nuc = 1 - C_P/(k_l_nuc + C_P);
        polymer_inh_growth = 1 - C_P/(k_l_growth + C_P);

        G_b = k_g *(S_b-1)^g * polymer_inh_growth;

        B_s_eff = B_s * polymer_inh_nuc;

        G_s = k_g *(S_s-1)^g * polymer_inh_growth;

        % Scaling of moment differential equations of surface

        dmu_s0dt_scaled = B_s_eff / mu_s0_char; % mu_s0 (t = 0) = 0
        dmu_s1dt_scaled = G_s * mu_s0 / mu_s1_char; % mu_s1 (t = 0) = 0
        dmu_s2dt_scaled = 2 * G_s * mu_s1 / mu_s2_char; % mu_s2 (t = 0) = 0
        dmu_s3dt_scaled = 3 * G_s * mu_s2 / mu_s3_char; % mu_s3 (t = 0) = 0

        % Calculating V_s

        %r_ddis = zhat*w_ccini^(1/3)*(w_cc)^(2/3)*(C_ds - C_db)/r_ccini;

        r_ddis = real(zhat*w_ccini^(1/3)*(w_cc)^(2/3)*(C_ds - C_db)/r_ccini);

        Vsingle_cc = w_cc/(rho*no_cc);
        r_cc = (Vsingle_cc*3/(4*pi))^(1/3);

        h_d = r_cc;

        A_p = no_cc * 4 * pi * r_cc ^ 2;
        V_s = A_p * h_d;

        % Scaling differential equation of w_cc
        % k_v, f_dsf, f_dcc are dimensionless; hence no scaling required

        dw_ccdt_scaled = (r_ddis + 3*rho*k_v*G_s*mu_s2*V_s*f_dsf)/(-f_dcc*w_cc_char);

        dw_ccdt = (r_ddis + 3*rho*k_v*G_s* mu_s2 *V_s*f_dsf)/(-f_dcc);

        % dV_s = d(A_p)*h_d+A_p*d(h_d)
        % d(A_p) = no_cc*8*pi*r_cc*d(r_cc)
        % w_cc = no_cc*4/3 * rho* pi * r_cc^3
        % d(w_cc) = no_cc*4*pi*rho*r_cc^2 *d(r_cc)
        % d(r_cc) = 1/(no_cc*4*pi*rho*r_cc^2) * d(w_cc)

        dV_sdt = no_cc*8*pi*r_cc*h_d* dw_ccdt/(no_cc*4*pi*rho*r_cc^2) + A_p* dw_ccdt/(no_cc*4*pi*rho*r_cc^2);

        if (h_d > 0.0030)
            h_d = 0.0030; % (cm)

            A_p = no_cc * 4 * pi * r_cc ^ 2;
            V_s = A_p * h_d;

            dV_sdt = no_cc*8*pi*r_cc*h_d* dw_ccdt/(no_cc*4*pi*rho*r_cc^2);
        end
        
        % Scaling of moment differential equations of bulk

        dmu_b0dt_scaled = B_b / mu_b0_char - (mu_s0 / mu_b0_char - mu_b0_scaled)*dV_sdt/(V_t - V_s); % mu_b0 (t = 0) = 0
        dmu_b1dt_scaled = G_b * mu_b0 / mu_b1_char - (mu_s1 / mu_b1_char - mu_b1_scaled)*dV_sdt/(V_t - V_s); % mu_b1 (t = 0) = 0
        dmu_b2dt_scaled = 2 * G_b * mu_b1 / mu_b2_char - (mu_s2 / mu_b2_char - mu_b2_scaled)*dV_sdt/(V_t - V_s); % mu_b2 (t = 0) = 0
        dmu_b3dt_scaled = 3 * G_b * mu_b2 / mu_b3_char - (mu_s3 / mu_b3_char - mu_b3_scaled)*dV_sdt/(V_t - V_s); % mu_b3 (t = 0) = 0

        dC_dbdt_scaled = (r_ddis/(V_t - V_s) - (C_ds - C_db_scaled * C_db_char)*dV_sdt/(V_t-V_s) - 3*rho*k_v*G_b*mu_b2*f_dsf) / C_db_char;

        dXdt_scaled = [dmu_s0dt_scaled; dmu_s1dt_scaled; dmu_s2dt_scaled; dmu_s3dt_scaled; dw_ccdt_scaled; dmu_b0dt_scaled; dmu_b1dt_scaled; dmu_b2dt_scaled; dmu_b3dt_scaled; dC_dbdt_scaled];

    end


end