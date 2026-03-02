function stats = Error_Analysis(p_exact, p_computed) 

% compute residual
residual = p_exact - p_computed;

% Compute Statistics
stats.infinityNorm = norm(residual, Inf); %infinity norm 
stats.mean_err = mean(residual); %mean error 
stats.variance = var(residual); % variance
stats.rms_err = sqrt(mean(residual.^2)); %root mean square error 
end