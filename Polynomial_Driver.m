%% Test Driver Summary
% -
%% Global Parameters
clear; clc; close all; 

% define functions
% f = @(x) (x - 2).^9 % function 1:  
% f =                 % function 2: product function parameterized by d
% f =                 % function 3: Lagrange Basis function parameterized by n
f = @(x) 1./ (1 + 25 * x.^2); % function 4: Runge's function 


% Define Interval [a,b] for Barycentric-2 coefficients 
Interval = [-1,1];
a = Interval(1);
b = Interval(2);

% Generate sample points for interpolation
numPoints = 100; % Number of points in the interval
xx = linspace(a, b, numPoints);
y = f(xx); % Evaluate the function at the sample points

% Define Precision ("single" or "double")
precision = "single";

% Polynomial Degree
n_list = [5, 10, 25, 50, 100];

% Mesh Type ("uniform", "chebychev1", or "chebychev2")
mesh_type = "uniform";

% Ordering ("increasing", "decreasing", or "leja")
order_type = "increasing"; 

n = 25;  % Test degree
k = 0:n; % Index from 0,...,n

% Generate mesh points from [-1,1] aka nodes for Interpolation
if mesh_type == "uniform"
    x = -1 + (2*k) / n;
elseif mesh_type == "chebychev1"
    x = cos((2*k + 1) * pi / (2*n +2)); 
elseif mesh_type == "chebychev2"
    x = cos(k * pi / n); 
end 

% Create mapping from [-1,1] -> [a,b]
mapped_mesh = 0.5*(a+b) + 0.5*(b-a) * x;  
ordered_mesh = Ordering(mapped_mesh, order_type); 

fprintf("Mesh: %s | Order: %s | n: %d\n", mesh_type, order_type, n);

%% Part 2
%% Part 3
%% Part 4
%% Part 5

