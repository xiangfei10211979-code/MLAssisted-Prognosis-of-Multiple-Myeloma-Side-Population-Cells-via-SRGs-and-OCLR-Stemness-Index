% Filename: Cox_Opt.m
function b_value = Cox_Opt(Current_TargetData, Current_W, Current_T, Current_censor)
    % Ensure W is a column vector for matrix multiplication
    if size(Current_W, 1) == 1 % If W is a row vector
        Current_W = Current_W'; % Convert to column vector
    end
    
    WeightedData = Current_TargetData * Current_W;
    
    % If the weighted data is constant, coxphfit may fail or return meaningless results
    if all(abs(diff(WeightedData)) < eps) || all(WeightedData == WeightedData(1))
        b_value = -realmax; % Return a very poor value because we aim to maximize b (minimize -b)
        return;
    end

    try
        [b_temp, ~, ~, ~] = coxphfit(WeightedData, Current_T, 'Censoring', Current_censor);
        if isempty(b_temp) || any(isnan(b_temp))
            % If b is empty or NaN, coxphfit failed
            b_value = -realmax; 
        else
            b_value = b_temp;
        end
    catch ME
        % If coxphfit throws an error
        fprintf('coxphfit error in Cox_Opt: %s\n', ME.message);
        b_value = -realmax; % Return a poor value when error occurs
    end
end
