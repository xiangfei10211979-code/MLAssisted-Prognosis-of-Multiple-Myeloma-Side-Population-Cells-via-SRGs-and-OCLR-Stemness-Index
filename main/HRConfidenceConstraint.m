% Filename: hrConfidenceIntervalConstraint.m
function [c, ceq] = hrConfidenceIntervalConstraint(W, TargetData_constr, T_constr, censor_constr)
    % Ensure W is a column vector
    if size(W, 1) == 1 % If W is a row vector
        W = W'; % Convert to column vector
    end
    
    NewData_constr = TargetData_constr * W;

    % If the weighted data is constant
    if all(abs(diff(NewData_constr)) < eps) || all(NewData_constr == NewData_constr(1))
        c = 1e6; % Return a very large positive number indicating strong constraint violation
        ceq = [];
        return;
    end

    try
        [b, ~, ~, stats] = coxphfit(NewData_constr, T_constr, 'Censoring', censor_constr);

        % Check if b and stats.se are valid
        if isempty(b) || any(isnan(b)) || isempty(stats.se) || any(isnan(stats.se)) || stats.se <= 0
            c = 1e6; % Strong constraint violation (e.g., invalid standard error)
        else
            z_score = 1.96; % Corresponds to 95% confidence interval
            
            b_lower = b - z_score * stats.se;
            b_upper = b + z_score * stats.se;

            hr_lower = exp(b_lower);
            hr_upper = exp(b_upper);
            
            % Check again whether hr_upper and hr_lower are valid (e.g., avoid Inf after exp())
            if isinf(hr_upper) || isnan(hr_upper) || isinf(hr_lower) || isnan(hr_lower)
                c = 1e6; % Invalid HR confidence interval, considered constraint violation
            else
                % Constraint c(W) <= 0  =>  (hr_upper - hr_lower) - 100 <= 0
                c = (hr_upper - hr_lower) - 100; 
            end
        end
    catch ME
        % Error handling if coxphfit fails
        fprintf('coxphfit error in hrConfidenceIntervalConstraint: %s\n', ME.message);
        c = 1e6; % Treat error as strong constraint violation
    end
    
    ceq = []; % No nonlinear equality constraints
end
