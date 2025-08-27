% ========================================================================
% Nonlinear Cantilever Beam Analysis using Ritz Method + Newton-Raphson
% ========================================================================
% Author   : Amal
% Date     : 27-08-2025
%
% Description:
%   This script performs geometrically nonlinear static analysis of a 
%   cantilever beam subjected to a tip load, using the Ritz method with 
%   polynomial trial functions. An explicit Newton-Raphson iterative solver 
%   is implemented for equilibrium iterations.
%
% Approximation:
%   - Axial displacement (u): 5th-order polynomial
%   - Transverse displacement (w): 3rd-order polynomial
%   - Single-element formulation (global Ritz approach)
%
% Outputs:
%   - Load vs. tip deflection curve
%   - Final deformed shape (linear vs. nonlinear)
% ========================================================================

clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Define Symbolic Variables and Functions
% -------------------------------------------------------------------------
syms x L E A I P real;                     % Beam/material parameters
syms c2 c3 d1 d2 d3 d4 d5 real;            % Ritz coefficients

% Trial functions for transverse (w) and axial (u) displacements
w_sym = c2*x^2 + c3*x^3;                   % Transverse displacement
u_sym = d1*x + d2*x^2 + d3*x^3 + d4*x^4 + d5*x^5; % Axial displacement

% Strains (von Kármán type nonlinear kinematics)
dw_dx_sym   = diff(w_sym, x);
d2w_dx2_sym = diff(w_sym, x, 2);
du_dx_sym   = diff(u_sym, x);

eps0_sym = du_dx_sym + 0.5 * dw_dx_sym^2;  % Axial strain
eps1_sym = -d2w_dx2_sym;                   % Bending curvature

% Strain energy density and total strain energy
U_density = 0.5 * (E*A * eps0_sym^2 + E*I * eps1_sym^2);
U_total   = int(U_density, x, 0, L);

%% ------------------------------------------------------------------------
% 2. Generate Residual Vector and Tangent Stiffness Matrix
% -------------------------------------------------------------------------
coeffs = [c2, c3, d1, d2, d3, d4, d5];     % Ritz unknowns
Pi_total = U_total;                        % Potential energy (internal)

fprintf('Generating symbolic residuals and tangent stiffness...\n');
R_int_expr = jacobian(Pi_total, coeffs)';  % Internal force vector
K_T_expr   = jacobian(R_int_expr, coeffs); % Tangent stiffness matrix

% Convert to fast MATLAB function handles
vars = [{L, E, A, I}, coeffs];
R_int_handle = matlabFunction(R_int_expr, 'Vars', vars);
K_T_handle   = matlabFunction(K_T_expr, 'Vars', vars);

%% ------------------------------------------------------------------------
% 3. Numerical Setup
% -------------------------------------------------------------------------
L_val = 1.0; E_val = 210e9; b_val = 0.05; h_val = 0.05;
A_val = b_val * h_val; 
I_val = b_val * h_val^3 / 12;

% Loading parameters
P_max = 5e5;                % Maximum tip load [N]
num_steps = 20;             % Number of load increments
load_steps = linspace(P_max/num_steps, P_max, num_steps);

% Newton-Raphson parameters
max_iter = 20;
tolerance = 1e-6;

%% ------------------------------------------------------------------------
% 4. Incremental-Iterative Newton-Raphson Procedure
% -------------------------------------------------------------------------
C = zeros(length(coeffs), 1);              % Initial Ritz coefficients
C_history = zeros(length(coeffs), num_steps);

for step = 1:num_steps
    current_load = load_steps(step);
    fprintf('\nLoad Step %d/%d, Applied Load = %.2e N\n', ...
        step, num_steps, current_load);
    
    % External force vector contribution (only w(L) terms)
    % P*w(L) -> derivative wrt c2 = P*L^2, wrt c3 = P*L^3
    F_ext = zeros(length(coeffs), 1);
    F_ext(1) = current_load * L_val^2;
    F_ext(2) = current_load * L_val^3;
    
    % Newton iterations
    for iter = 1:max_iter
        % Current coefficient values
        c_vals = num2cell(C');
        
        % Internal forces and tangent stiffness
        F_int = R_int_handle(L_val, E_val, A_val, I_val, c_vals{:});
        K_T   = K_T_handle(L_val, E_val, A_val, I_val, c_vals{:});
        
        % Residual = Internal - External
        R = F_int - F_ext;
        norm_R = norm(R);
        fprintf('  Iteration %d: Residual Norm = %.4e\n', iter, norm_R);
        
        % Convergence check
        if norm_R < tolerance
            fprintf('  Converged in %d iterations.\n', iter);
            break;
        end
        if iter == max_iter, warning('Max iterations reached.'); end
        
        % Update coefficients
        delta_C = -K_T \ R;
        C = C + delta_C;
    end
    
    C_history(:, step) = C; % Save step history
end

%% ------------------------------------------------------------------------
% 5. Post-Processing: Load-Deflection Curve
% -------------------------------------------------------------------------
tip_deflections = zeros(1, num_steps);
for i = 1:num_steps
    c2_val = C_history(1, i);
    c3_val = C_history(2, i);
    tip_deflections(i) = c2_val * L_val^2 + c3_val * L_val^3;
end

figure('Name', 'Load-Deflection Curve');
plot(tip_deflections, load_steps, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Nonlinear Ritz-Newton');
hold on;

% Linear elastic closed-form deflection
w_linear = (load_steps * L_val^3) / (3 * E_val * I_val);
plot(w_linear, load_steps, '--r', 'LineWidth', 1.5, ...
    'DisplayName', 'Linear Theory');

grid on; xlabel('Tip Deflection w (m)'); ylabel('Load P (N)');
title('Load-Deflection Curve (Cantilever Beam)');
legend('show', 'Location', 'NorthWest');
set(gca, 'FontSize', 12);

%% ------------------------------------------------------------------------
% 6. Final Deformed Shape
% -------------------------------------------------------------------------
x_points = linspace(0, L_val, 200);

% Nonlinear solution (final step coefficients)
w_sol_final(x) = subs(w_sym, coeffs, C');
u_sol_final(x) = subs(u_sym, coeffs, C');
w_nonlinear = double(w_sol_final(x_points));
u_nonlinear = double(u_sol_final(x_points));
x_deformed  = x_points + u_nonlinear;

% Linear solution
P_final = current_load;
w_linear = (P_final * x_points.^2) ./ (6 * E_val * I_val) .* (3*L_val - x_points);

% Plot deformed shapes
figure('Name', 'Final Deformed Shape');
plot([0, L_val], [0, 0], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Undeformed');
plot(x_points, w_linear, 'r:', 'LineWidth', 2, 'DisplayName', 'Linear');
plot(x_deformed, w_nonlinear, 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear');
grid on; axis equal;
xlabel('Axial Position x (m)'); ylabel('Transverse Deflection w (m)');
title(sprintf('Final Deformed Shape at P = %.2e N', P_final));
legend('show', 'Location', 'NorthWest');
set(gca, 'FontSize', 12);
