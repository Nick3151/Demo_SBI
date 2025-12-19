function metrics = trial_balance_metrics(n_correct, n_miss, n_error)

    % Total number of trials
    total = n_correct + n_miss + n_error;

    % Proportions
    props = [n_correct, n_miss, n_error] / total;

    % Entropy-based balance (base 3 for normalization between 0 and 1)
    entropy_balance = -sum(props .* log(props) / log(3), 'omitnan');

    % Skew index: max deviation from uniform distribution (1/3)
    skew_index = max(abs(props - 1/3));

    % Gini balance: 1 - sum of squared proportions
    gini_balance = 1 - sum(props .^ 2);

    % Return as a struct
    metrics.entropy_balance = round(entropy_balance, 3);
    metrics.skew_index = round(skew_index, 3);
    metrics.gini_balance = round(gini_balance, 3);

end
