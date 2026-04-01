%% Test Driver Summary
% -
%% Global Parameters 
clear; clc; close all; 

% Options # 1-4 
function_choice = 4; 

% Define Precision ("single" or "double")
precision = "single";

% Mesh Type ("uniform", "chebychev1", or "chebychev2")
mesh_type = "uniform";

% Ordering ("increasing", "decreasing", or "leja")
order_type = "increasing"; 

% Method ("newton", "barycentric1", or "barycentric2")
method = "barycentric2"; 


n = 10; 
f = @(x) EvaluateFunction(x, function_choice, n); 

%% Mesh Generation 

%Generate Mesh
k = 0:n; % Index from 0,...,n

% Generate mesh points from [-1,1] aka nodes for Interpolation
if mesh_type == "uniform"
    x = -1 + (2*k) / n;
elseif mesh_type == "chebychev1"
    x = cos((2*k + 1) * pi / (2*n +2)); 
elseif mesh_type == "chebychev2"
    x = cos(k * pi / n); 
end 

% Define Interval [a,b] for Barycentric-2 coefficients 
Interval = [-1,1];
a = Interval(1);
b = Interval(2);

% Create mapping from [-1,1] -> [a,b]
mapped_mesh = 0.5*(a+b) + 0.5*(b-a) * x;  

% Generate ordered mesh based on the specified ordering type
ordered_mesh = Ordering(mapped_mesh, order_type); 

% Use product form as exact reference for f1 and f2, double precision otherwise
if function_choice == 3
    if method == "barycentric2"
        % Build y_double against the xi mesh that Barycentric2 will use internally
        [~, ~, xi_temp] = Barycentric2Coefficients(Interval, f, mesh_type, n, "double");
        y_double = zeros(size(xi_temp));
        y_double(end) = 1;
    else
        y_double = zeros(size(ordered_mesh));
        y_double(end) = 1;
    end

 % evaluated in product form to avoid cancellation
elseif function_choice == 1 || function_choice == 2
    y_double = EvaluateProductForm(ordered_mesh, function_choice, n);
else
    y_double = f(ordered_mesh);
end

% Generate sample points — must be defined BEFORE y is evaluated
numPoints = 100;
xx = linspace(a, b, numPoints);
xx = xx(:);

if function_choice == 1 || function_choice == 2
    y = EvaluateProductForm(xx, function_choice, n);
elseif function_choice == 3
    % Exact f3 must be evaluated as the actual Lagrange basis function
    % Use double precision Barycentric2 as the reference
    [beta_d, ~, xi_d] = Barycentric2Coefficients(Interval, f, mesh_type, n, "double");
    y_double_ref = zeros(size(xi_d)); y_double_ref(end) = 1;
    y = Barycentric2_Interpolation(xx, xi_d, beta_d, y_double_ref, "double");
else
    y = f(xx);
end


%% Conditioning functions calls
[px_ref, k1, k2] = Barycentric_Conditioning(xx, ordered_mesh, y_double);
Lambda_n = max(k1); 
H_n = max(k2); 

fprintf('--- Conditioning Summary (n=%d) ---\n', n);
fprintf('Lebesgue Constant (Lambda_n): %.4e\n', Lambda_n);
fprintf('Condition Summary (H_n):      %.4e\n\n', H_n);

fprintf("Mesh: %s | Order: %s | n: %d\n", mesh_type, order_type, n);

%% Coefficient and Interpolation function calls

% Cast node data to the requested precision, then compute interpolation
% coefficients and evaluate the interpolant at the sample points.

% Method Type ("barycentric1", "barycentric2", or "newton")
x_single = cast(ordered_mesh, 'single');
y_single = cast(y_double, 'single');

if method == "barycentric1"
    % Weights gamma_j = 1 / prod_{k != j}(x_j - x_k)
    % Interpolant: p(x) = w(x) * sum_j [ gamma_j * f_j / (x - x_j) ]
    %   where w(x) = prod_j(x - x_j)
    
    [gamma, fx] = Barycentric1Coefficients(ordered_mesh, f , precision);
    px = Barycentric1_Interpolation(xx, ordered_mesh, gamma, fx, precision); 

elseif method == "barycentric2"
    % Simplified weights beta_j with analytic formulas per mesh type.
    % Interpolant: p(x) = [ sum_j beta_j*f_j/(x-x_j) ] / [ sum_j beta_j/(x-x_j) ]

    [beta, fx, xi] = Barycentric2Coefficients(Interval, f, mesh_type, n, precision);
    if function_choice == 3
        %override function values with exact lagrange basis 
        fx = cast(y_double, precision);
    end
    px = Barycentric2_Interpolation(xx, xi, beta, fx, precision);

elseif method == "newton"
    % Coefficients are the leading diagonal of the divided-difference table.
    % Interpolant evaluated efficiently via Horner's method.

    ordered_mesh_cast = cast(ordered_mesh, precision);
    if function_choice == 3

        % Build divided differences directly from y_double node values
        % calling f(), which would return zeros everywhere.

        fx = cast(y_double, precision);
        newton_coefficients = double(fx);

        for j = 2:length(ordered_mesh_cast)
            for i = length(ordered_mesh_cast):-1:j
                % In-place divided difference update (O(n^2) algorithm)
                newton_coefficients(i) = (newton_coefficients(i) - newton_coefficients(i-1)) / (ordered_mesh_cast(i) - ordered_mesh_cast(i-j+1));
            end
        end

        % Standard path: evaluate f at nodes and compute divided differences
        newton_coefficients = cast(newton_coefficients, precision);
    else
        [fx, newton_coefficients] = NewtonCoefficients(ordered_mesh_cast, f, precision);
    end
    px = Newton_Interpolation(xx, ordered_mesh_cast, newton_coefficients, precision);
end

% Ensure outputs are column vectors for consistent downstream arithmetic
px = px(:); 
y = y(:);

fprintf('Coefficients computed using %s method\n', method);

%% Error Analysis Function Calls

% Compare the interpolant p_n(x) against the exact reference y(x) at all
% sample points, and also verify that p_n interpolates exactly at the nodes.

stats = Error_Analysis(y, double(px)); 

% Relative infinity-norm error (normalises by the scale of f)
rel_error = norm(y - double(px), Inf) / norm(y, Inf);

% Node interpolation check: p_n(x_i) should equal f(x_i) to machine precision.
% Only meaningful for Newton 
if method == "newton"
    node_vals = Newton_Interpolation(ordered_mesh_cast, ordered_mesh_cast, newton_coefficients, precision);
    node_error = norm(double(node_vals(:)) - double(fx(:)), Inf);
else
    node_error = NaN; % not applicable for barycentric methods
end

fprintf('\n--- Statistical Results (n = %d) ---\n', n); 
fprintf('Infinity Norm: %.4e\n', stats.infinityNorm);
fprintf('Relative Infinity Norm: %.4e\n', rel_error);
fprintf("RMS Error: %.4e\n", stats.rms_err); 
fprintf("Mean Error: %.4e\n", stats.mean_err); 
fprintf("Error Variance: %.4e\n", stats.variance); 
fprintf('Node interpolation error: %.4e\n', node_error);

% Diagnostic: confirm data types used in computation
fprintf('fx class: %s, size: %d\n', class(fx), length(fx));
if method == "newton"
    fprintf('newton_coefficients class: %s\n', class(newton_coefficients));
    fprintf('ordered_mesh_cast class: %s\n', class(ordered_mesh_cast));
end
fprintf('ordered_mesh class: %s\n', class(ordered_mesh));

%% Interpolation Visualizations 

% Figure 1: (top) Exact function, interpolant, and node markers overlaid;
% (bottom) Absolute residual |f(x) - p_n(x)| on a log scale.

figure("Name","Interpolation Analysis"); 
subplot(2,1,1); 
plot(xx, y, 'k-', 'LineWidth', 2); hold on;
plot(xx, px, 'r--', 'LineWidth', 1.5);
plot(ordered_mesh, fx, 'bo');
title(sprintf('f%d(x) | %s mesh | %s order | n=%d | method=%s', ...
    function_choice, mesh_type, order_type, n, method));
xlabel('x'); ylabel('p(x)');
legend('Exact f(x)', 'Interpolant p_n(x)', 'Nodes', 'Location', 'best');
grid on;

subplot(2,1,2);
semilogy(xx, abs(y(:) - px(:)), 'm', 'LineWidth', 1.5);
title(sprintf('Absolute Residual |f(x) - p_n(x)|  ||r||_inf = %.4e', stats.infinityNorm));
xlabel('x'); ylabel('|r(x)|');
grid on;


%% Conditioning Visualizations 

% Figure 2: (top) Lebesgue function kappa_1(x) and its maximum Lambda_n;
% (bottom) Condition number kappa_2(x) and its maximum H_n.


figure("Name", "Conditioning Analysis");

subplot(2,1,1);
semilogy(xx, k1, 'b-', 'LineWidth', 1.5); hold on;
semilogy(xx, ones(size(xx)) * Lambda_n, 'r--', 'LineWidth', 1.2);
title(sprintf('Lebesgue Function \\kappa(x,n,1) | %s mesh | n=%d', mesh_type, n));
xlabel('x'); ylabel('\kappa(x,n,1)');
legend('\kappa(x,n,1)', sprintf('\\Lambda_n = %.4e', Lambda_n), 'Location', 'best');
grid on;

subplot(2,1,2);
semilogy(xx, k2, 'r-', 'LineWidth', 1.5); hold on;
semilogy(xx, ones(size(xx)) * H_n, 'b--', 'LineWidth', 1.2);
title(sprintf('Condition Number \\kappa(x,n,y) | %s mesh | n=%d', mesh_type, n));
xlabel('x'); ylabel('\kappa(x,n,y)');
legend('\kappa(x,n,y)', sprintf('H_n = %.4e', H_n), 'Location', 'best');
grid on;


%% Task 5 — Convergence Testing for f4 only
% For the Runge function, sweep over multiple values of n and all three mesh
% types to visualise how the approximation error decays as n increases.
% Barycentric Form 2 is used throughout for consistency and stability.

if function_choice == 4
    n_values = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100];
    mesh_types = {"uniform", "chebychev1", "chebychev2"};
    colors = {'r-o', 'b-s', 'g-^'};
    all_errors = zeros(length(mesh_types), length(n_values));

    for m = 1:length(mesh_types)
        mt = mesh_types{m};
        for idx = 1:length(n_values)
            ni = n_values(idx);
            f4 = @(x) EvaluateFunction(x, 4, ni);

            % Build the Barycentric2 interpolant for this (mesh type, degree)
            [beta_i, fi_b, xi] = Barycentric2Coefficients([-1,1], f4, mt, ni, 'single');
            xi = xi(:);

            % Evaluate interpolant on a grid for error measurement
            xx_fine = linspace(-1, 1, 300)';
            px_fine = Barycentric2_Interpolation(xx_fine, xi, beta_i, fi_b, 'single');
            y_fine = EvaluateFunction(xx_fine, 4, ni);

            % Record the infinity-norm error for this (mesh type, degree) pair
            all_errors(m, idx) = norm(double(px_fine(:)) - y_fine(:), Inf);
        end
    end

    % Print table for convergence summary 
    fprintf('\n--- Convergence Sweep (f4, Barycentric Form 2) ---\n');
    fprintf('%-6s %-12s %-12s %-12s\n', 'n', 'Uniform', 'Cheb1', 'Cheb2');
    for idx = 1:length(n_values)
        fprintf('%-6d %-12.4e %-12.4e %-12.4e\n', n_values(idx), ...
            all_errors(1,idx), all_errors(2,idx), all_errors(3,idx));
    end

    % Figure 3: Convergence curves for each mesh type
    figure("Name", "Convergence Sweep — f4(x)");
    hold on;
    for m = 1:length(mesh_types)
        semilogy(n_values, all_errors(m,:), colors{m}, 'LineWidth', 1.5, 'MarkerSize', 6);
    end
    yline(1e-7, 'k--', 'LineWidth', 1.2);
    title('Convergence of p_n(x) \rightarrow f_4(x) as n increases');
    xlabel('n'); ylabel('||f - p_n||_\infty');
    legend('Uniform', 'Chebyshev 1st kind', 'Chebyshev 2nd kind', ...
           'Single precision floor', 'Location', 'best');
    grid on;
    set(gca, 'XTick', n_values);
end


% =========================================================================
%                         FUNCTIONS
%=========================================================================
%% Barycentric 1 Coefficients function

% Computes the Barycentric Form 1 weights gamma_j for a given mesh.

function [gamma, fx] = Barycentric1Coefficients(mesh, f, precision)

    n = length(mesh) - 1;
    mesh = double(mesh); % Compute weights in high precision
    gamma = ones(1, n+1); % Store Coefficients/weights

    % Build gamma_j = prod_{k != j} (x_j - x_k) for each node j
    for j = 1:n+1
        for k = 1:n+1
            if j ~= k
                gamma(j) = gamma(j) * (mesh(j) - mesh(k));
            end
        end
    end
    gamma = 1 ./ gamma; % Convert to weights 

    % Defined by single and double point precision 
    fx = cast(f(mesh), precision);
    gamma = cast(gamma, precision);
end



%% Barycentric 1 Interpolation 
function px = Barycentric1_Interpolation(evaluation_points, mesh, gamma, fx, precision)

% Initialize input values as single or double floating point precision
if precision == "single"
    evaluation_points = single(evaluation_points);
    mesh = single(mesh);
    gamma = single(gamma); 
    fx = single(fx);
else
    mesh = double(mesh);
    evaluation_points = double(evaluation_points);
    gamma = double(gamma); 
    fx = double(fx); 
end

% Initialize output array with the appropriate type
px = cast(zeros(size(evaluation_points)), precision);

% loop through each point that gets evaluated 
for j=1:length(evaluation_points)
    x_current = evaluation_points(j); 

    % Compute difference vector (x - x_i)
    difference = x_current - mesh; 

    % check if (x-x_i) = 0 at any point
    % if difference = 0  => Code fails via (0*inf) 
    exact_spot = find(difference == 0,1);

    % if exact_spot is Not empty => return f(x) 
    if ~isempty(exact_spot)
        px(j) = fx(exact_spot);
    % Otherwise compute Barycentric 1 formula
    else 
        lx = prod(difference); 
        summation = sum((gamma .* fx ./ difference));
        px(j) = lx * summation; 
    end
end
end


%% Barycentric 2 Coefficients
function [beta, fx, x_mesh] = Barycentric2Coefficients(Interval, f, mesh_type, n, precision)
a = Interval(1); 
b = Interval(2); 

% Index from 0,...,n 
k = 0:n; 

% Generate mesh points from [-1,1] aka nodes for Interpolation
if mesh_type == "uniform"
    x = -1 + (2*k) / n;
elseif mesh_type == "chebychev1"
    x = cos((2*k + 1) * pi / (2*n +2)); 
elseif mesh_type == "chebychev2"
    x = cos(k * pi / n); 
end 

% Create mapping from [-1,1] -> [a,b]
x_mesh = 0.5*(a+b) + 0.5*(b-a) * x;  

% Initialize precision for mesh and function eval.
x_mesh = cast(x_mesh, precision); 
fx = cast(f(x_mesh), precision);

% Compute Coefficients for Barycentric form 2 
if mesh_type == "uniform"
    beta = cast(zeros(1,n+1), precision);
    for j = 0:n
        % Compute binomial(n,j) in log space to avoid integer overflow
        beta(j+1) = (-1)^j * exp(gammaln(n+1) - gammaln(j+1) - gammaln(n-j+1));
    end

% if chebychev points of 1st kind => Sinusodial Coefficients/weights 
elseif mesh_type == "chebychev1"
    beta = ((-1).^k) .* sin((2*k + 1) * pi / (2*n + 2));

% if chebychev points of 2nd kind => [1/2, 1, 1, ..., 1, 1, 1/2] * (-1)^j   
elseif mesh_type == "chebychev2"
    beta = cast(ones(1, n+1), precision);
    beta(1) = 0.5; 
    beta(end) = 0.5;
    beta = beta .* ((-1).^k);
end 

% Ensure beta is computed with correct precision
beta = cast(beta, precision); 
end



%% Barycentric 2 Interpolation 
function px = Barycentric2_Interpolation(evaluation_points, x_mesh, beta, fx, precision)

% Initialize Sroage for polynomial construction
px = cast(zeros(size(evaluation_points)), precision);

%Cast Single/double point precision on input parameters
evaluation_points = cast(evaluation_points(:), precision);  
x_mesh = cast(x_mesh(:)', precision);   
beta = cast(beta(:)', precision);       
fx = cast(fx(:)', precision);           

% Begin polynomial evaluation
for j=1:length(evaluation_points)
    x_current = evaluation_points(j); 

    % Compute difference vector (x - x_i) with transformed mesh
    difference = x_current - x_mesh; 

    % check if (x-x_i) = 0 at any point
    % if difference = 0  => Code fails via (0*inf) 
    exact_spot = find(difference == 0,1);

    % if exact_spot is Not empty => return f(x) 
    if ~isempty(exact_spot)
        px(j) = fx(exact_spot);
    % Otherwise compute Barycentric 2 formula
    else
        rho = beta ./ difference;
        px(j) = sum(rho .* fx)/ sum(rho);
    end 
end
end



%% Divided Difference Function - O(n^2) 
% function [fx, DDT, newton_coefficients] = Divided_Difference(mesh, f, precision)
% 
% % Num of nodes
% n = length(mesh); 
% % Evaluate the function at the mesh points
% fx = cast(f(mesh), precision);   
% % Storage for Table with same precision as fx 
% DDT = zeros(n, n, 'like', fx); 
% 
% for i = 1:n 
%     DDT(1,i) = fx(i); 
% end
% 
% for j = 2:n
%     for i = 1:(n - j + 1)
%         num = DDT(j-1, i+1) - DDT(j-1, i);
%         den = mesh(i + j - 1) - mesh(i);
%         DDT(j, i) = num / cast(den, precision);
%     end
% end
% 
% newton_coefficients = DDT(:,1); 
% end



%% Newton Coefficients Function - O(n)

function [fx, newton_coefficients] = NewtonCoefficients(mesh, f, precision)
    n = length(mesh); 
    fx = cast(f(mesh), precision);   
    newton_coefficients = double(fx); % Accumulate in double to reduce rounding

    % In-place update: after the j-th outer iteration, newton_coefficients(j:end)
    % holds divided differences of order (j-1) anchored at successive nodes.
    for j = 2:n
        for i = n:-1:j
            newton_coefficients(i) = (newton_coefficients(i) - newton_coefficients(i-1)) / (mesh(i) - mesh(i-j+1));
        end
    end
    newton_coefficients = cast(newton_coefficients, precision);
end



%% Newton Interpolation Function 

function px = Newton_Interpolation(evaluation_points, mesh, newton_coefficients, precision)

    % Cast all inputs to the working precision
evaluation_points = cast(evaluation_points, precision);
mesh = cast(mesh, precision);
coefficients  = cast(newton_coefficients, precision);

n = length(coefficients) - 1; 
px = cast(zeros(size(evaluation_points)), precision);

% Loop through each point we want to evaluate
for j = 1:length(evaluation_points)
    x = evaluation_points(j);
        
    % Horner's Method: Start from the highest order coefficient
    s = coefficients(n + 1); 
        
    % Iterate backwards from n-1 down to 0
    % s = s*(x - x_i) + alpha_i
    for i = n:-1:1
        s = s * (x - mesh(i)) + coefficients(i);
    end
        
    px(j) = s;
end
end



%% Ordering Function 

% Reorders the interpolation nodes according to the specified strategy.
%
%   "increasing" — ascending sort
%   "decreasing" — descending sort
%   "leja"       — Leja ordering: greedily maximises the product
%                  prod_{k < j} |x_j - x_k| at each step, which improves
%                  the numerical stability of Newton interpolation
%

function ordered_mesh = Ordering(mesh, type)
    n = length(mesh);
    if type == "increasing"
        ordered_mesh = sort(mesh, "ascend");
    elseif type == "decreasing"
        ordered_mesh = sort(mesh, "descend");
    elseif type == "leja"
        ordered_mesh = zeros(size(mesh));
        remaining = mesh;

        % Start with the node of largest absolute value
        [~, idx] = max(abs(remaining));
        ordered_mesh(1) = remaining(idx);
        remaining(idx) = [];

         % select the remaining nodes and at each step pick the node
        % that maximises product 
        for k = 2:n
            [~, idx] = max(prod(abs(remaining' - ordered_mesh(1:k-1)), 2));
            ordered_mesh(k) = remaining(idx);
            remaining(idx) = [];
        end
    else
        % Unknown type — return mesh unchanged
        ordered_mesh = mesh;
    end
end



%% Error Analysis Function
% Computes summary statistics for the pointwise error between exact and
% computed values.

function stats = Error_Analysis(p_exact, p_computed) 

% compute residual
residual = p_exact - p_computed;

% Compute Statistics
stats.infinityNorm = norm(residual, Inf); %infinity norm 
stats.mean_err = mean(residual); %mean error 
stats.variance = var(residual); % variance
stats.rms_err = sqrt(mean(residual.^2)); %root mean square error 
end



%% Conditioning Function 

% Computes the interpolant and two pointwise conditioning measures at each
% evaluation point using the Barycentric Form 1 representation.

function [px, kappa1, kappa2] = Barycentric_Conditioning(evaluation_points, mesh, fx)
    % Work in double precision for accurate conditioning estimates
    xx = double(evaluation_points(:));   
    x_i = double(mesh(:));               
    f_i = double(fx(:));                 
    n = length(x_i)-1; 

    % stoage for weights 
    gamma = ones(1, n+1); 

    % Compute Barycentric Form 1 weights gamma_j = 1/prod_{k!=j}(x_j - x_k)
    for j = 1:n+1 
        for k = 1:n+1 
            if j ~= k
                gamma(j) = gamma(j) * (x_i(j)-x_i(k)); 
            end
        end
    end
    gamma = 1 ./ gamma;
    gamma = gamma(:);                    

    % Pre-allocate outputs
    px = zeros(size(xx)); 
    kappa1 = zeros(size(xx));
    kappa2 = zeros(size(xx));

    for j = 1:length(xx)
        x = xx(j);
        difference = x - x_i;           

        % Handle the case where x coincides with a node
        exact_spot = find(difference == 0, 1);
        if ~isempty(exact_spot)
            % Interpolant equals the node value; conditioning is trivially 1
            px(j) = f_i(exact_spot); 
            kappa1(j) = 1.0; 
            kappa2(j) = 1.0; 
        else 
            % w(x) = prod_j (x - x_j)
            lx = prod(difference); 
            coefficients = gamma ./ difference;   % gamma_j / (x - x_j) for each j

            % Interpolant: p(x) = w(x) * sum(gamma_j * f_j / (x - x_j))
            px(j) = lx * sum(coefficients .* f_i);

            % Lebesgue function: sum of |gamma_j / (x - x_j)| weighted by |w(x)|
            kappa1(j) = abs(lx) * sum(abs(coefficients)); 

            % Relative condition number: |w(x)| * sum|coeff*f| / |p(x)|
            numerator = abs(lx) * sum(abs(coefficients .* f_i)); 
            kappa2(j) = numerator / abs(px(j)); 
        end
    end
end



%% Product form evaluation function

% Evaluates f1 and f2 via explicit multiplication rather than using the
% expanded polynomial form, which suffers from catastrophic cancellation
% near the clustered roots.

function y = EvaluateProductForm(x, function_choice, n)
    x = x(:);
    switch function_choice
        case 1
            % f1(x) = (x-2)^9 — single root, product form
            y = ones(size(x));
            for i = 1:9
                y = y .* (x - 2);
            end

        case 2
            % f2(x;n) = prod(x-i) for i=1..n
            y = ones(size(x));
            for i = 1:n
                y = y .* (x - i);
            end

        otherwise
            % Fall back to double precision for f3, f4
            y = double(EvaluateFunction(x, function_choice, n));
    end
end

%% Function Evaluation function 

function [y] = EvaluateFunction(x, type, n)
    switch type
        case 1 % function 1
            x = x(:);
            y = (x - 2).^9;

        case 2 % function 2
            x = x(:);
            y = ones(size(x));
                for i = 1:n
                    y = y .* (x - i);
                end

         case 3 % Lagrange basis — f(x) used only for sample plot, actual y_double set in driver
            x = x(:);
            y = zeros(size(x)); % returns 0 everywhere; driver overrides y_double at nodes

        case 4 % Runge function
            x = x(:);
            y = 1.0 ./ (1.0 + 25.0 .* x.^2); 
    end
end







