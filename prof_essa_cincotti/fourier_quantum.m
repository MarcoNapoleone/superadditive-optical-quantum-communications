% Clear all variables from workspace
clear all;

% Define a range of SNR values (nR) in logarithmic scale
nR = logspace(-5, 0, 100);  % 100 points between 10^-5 and 10^0
MM = 3;  % Parameter defining modulation order
M = 2^MM;  % Total number of states (modulation level)

% Auto-correlation term calculation
AUTO = nR * (2 * M^2 + 1) / (3 * M);

% Exponential decay for the auto-correlation term
INTENSITYTOTAL = exp(-AUTO);

% Entropy of the auto-correlation exponential decay
ENTROPYTOTAL = entropy(exp(-AUTO));

% Loop to accumulate cross-correlation contributions
for k = 1:M-1
    % Cross-correlation term
    CROSS = nR / (M * (sin(pi * k / M))^2);
    
    % Exponential decay for cross-correlation term
    CROSSEXPONENTIAL = exp(-CROSS);
    
    % Entropy for cross-correlation exponential decay
    CROSSENTROPY = entropy(CROSSEXPONENTIAL);
    
    % Add cross-correlation intensity and entropy to totals
    INTENSITYTOTAL = INTENSITYTOTAL + CROSSEXPONENTIAL;
    ENTROPYTOTAL = ENTROPYTOTAL + CROSSENTROPY;
end

% Calculate auto-correlation and cross-correlation for MM-modulation
AUTO1 = nR * (2 * MM^2 + 1) / (3 * MM^2);
CROSS1 = nR / (MM^2 * (sin(pi / MM))^2);

% Exponential terms for CROSS1 and AUTO1
CROSS1EXPONENTIAL = (1 + 2 * exp(-CROSS1) + exp(-AUTO1) + exp(-2 * CROSS1) ...
    + 2 * exp(-AUTO1 - CROSS1) + exp(-AUTO1 - 2 * CROSS1)) / 8;

% Calculate entropy contributions for different terms
CROSS1ENTROPY = 2 * entropy(exp(-CROSS1)) + entropy(exp(-AUTO1)) + ...
    entropy(exp(-2 * CROSS1)) + 2 * entropy(exp(-AUTO1 - CROSS1)) + ...
    entropy(exp(-AUTO1 - 2 * CROSS1));

% PF Modulation calculation
PFM = entropy(CROSS1EXPONENTIAL) - CROSS1ENTROPY / 8;

% Calculate different performance metrics
Holevo = entropy(1/2 - exp(-2 * nR) / 2);   % Holevo bound
Homodyne = 1 - entropy((1 + erf(sqrt(2 * nR))) / 2);  % Homodyne detection
Helstrom = 1 - entropy((1 - sqrt(1 - exp(-4 * nR))) / 2);  % Helstrom bound
GM = log(M) * (1 - exp(-M * nR)) / M;  % Gaussian Modulation (GM)
FM = entropy(INTENSITYTOTAL / M) - ENTROPYTOTAL / M;  % Final modulation (FM)

% Plot performance metrics on a log-log scale

% Figure 1: Comparison of different methods
figure(1)
loglog(nR, Holevo, '-k', 'DisplayName', 'Holevo');  % Plot Holevo bound
hold on
loglog(nR, Homodyne, '-g', 'DisplayName', 'Homodyne');  % Plot Homodyne detection
loglog(nR, Helstrom, ':g', 'DisplayName', 'Helstrom');  % Plot Helstrom bound
loglog(nR, GM, 'r', 'DisplayName', 'GM');  % Plot Gaussian Modulation
loglog(nR, FM, ':r', 'DisplayName', 'FM');  % Plot Final modulation
grid on
legend show;  % Show legend

% Figure 2: Focused comparison of Holevo, FM, and PFM
figure(2)
loglog(nR, Holevo, '-k', 'DisplayName', 'Holevo');  % Plot Holevo bound
hold on
loglog(nR, FM, ':r', 'DisplayName', 'FM');  % Plot Final modulation
loglog(nR, PFM, ':c', 'DisplayName', 'PFM');  % Plot PF Modulation
grid on
legend show;  % Show legend
