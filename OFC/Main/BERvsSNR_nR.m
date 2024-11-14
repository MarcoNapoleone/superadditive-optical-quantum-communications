clear;
close all;

%% Parameters
M = 32; % Fixed M value for Hadamard/Fourier Machine

% Array of mean photon numbers
n_R = [0.01,0.1, 1]

%% Channel parameters
SNRdB = 0:1:20; % SNR in dB array
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale

%% Preallocate arrays for BER
BER_teo_H = zeros(length(n_R), length(SNR)); % Theoretical BER for Hadamard (HM)
BER_teo_F = zeros(length(n_R), length(SNR)); % Theoretical BER for Fourier (FM)
BER_BPSK = zeros(length(M), length(SNR)); % Theoretical BER for BPSK

%% Loop over each n_R and SNR values
for nR_idx = 1:length(n_R)
    current_n_R = n_R(nR_idx); % Select current mean photon number

    for idx = 1:length(SNR)
        snr = SNR(idx);

        % Calculate noise photons for Green and Fourier machines
        n_N_H = current_n_R / snr; % Noise for Hadamard
        n_N_F = current_n_R / (snr * M); % Noise for Fourier

        %% Theoretical BER (Hadamard and Fourier Machines)
        BER_teo_H(nR_idx, idx) = BER_theoretical_Hadamard(current_n_R, n_N_H, M);
        BER_teo_F(nR_idx, idx) = BER_theoretical_Fourier(current_n_R, n_N_F, M);
        BER_BPSK(nR_idx, idx) = BER_the_BPSK(current_n_R, n_N_H, M);
    end

end

%% Plot Results for each n_R
f = figure;
f.Position = [10 10 550 550]; 
styles = {'-', '-.', ':'}; % Line styles for different M
colors = {'#098bf8', '#fe5f55', '#f7b801', '#2ab7ca', '#111111'}; % Colors for Hadamard, Fourier, BPSK

% Adjust axes position to leave space for the legend outside the plot
ax = gca;
ax.Position = [0.1, 0.2, 0.82, 0.63];  % Adjust axes position to give space above for the legend

% Loop through each n_R value
for nR_idx = 1:length(n_R)
    current_n_R = n_R(nR_idx);

    % Plot Hadamard results for current n_R
    semilogy(SNRdB, BER_teo_H(nR_idx, :), 'Color', colors{1}, 'LineStyle', styles{nR_idx}, ...
             'LineWidth', 1.5, 'DisplayName', 'Hadamard nR=' + string(current_n_R));
    hold on;

    % Plot Fourier results for current n_R
    semilogy(SNRdB, BER_teo_F(nR_idx, :), 'Color', colors{2}, 'LineStyle', styles{nR_idx}, ...
             'LineWidth', 1.5, 'DisplayName', 'Fourier nR=' + string(current_n_R));
    
    % Plot BPSK results for current n_R
    semilogy(SNRdB, BER_BPSK(nR_idx, :), 'Color', colors{3}, 'LineStyle', styles{nR_idx}, ...
             'LineWidth', 1.5, 'DisplayName', 'BPSK nR=' + string(current_n_R));
end

% Set axis labels and grid
xlabel('SNR (dB)');
ylabel('SER');
grid on;

% Adjust legend to be above the plot, with 3 columns, and adjust the size and position
lgd = legend('NumColumns', 3, 'Location', 'southoutside');
lgd.Position(2) = 0.85;  % Adjust the vertical position to ensure consistency

% Set consistent axes limits if needed
xlim([min(SNRdB), max(SNRdB)]);
ylim([1e-3, 1.5]);  % Adjust as needed

% Export the figure as a PDF
exportgraphics(gcf, 'output/BERvsSNR_nR.pdf');


%% Theoretical BER Calculation Functions

% BPSK theoretical BER calculation
function BER = BER_the_BPSK(n_R, n_N, M)
    BER = 1- (1 - (1-exp(-n_N)+exp(-n_R-n_N))/2) ^ log2(M);
end


% Hadamard (Green Machine) theoretical BER calculation
function BER = BER_theoretical_Hadamard(n_R, n_N, M)
    BER = (exp(- (M * n_R + n_N)) + ((M - 1) * (1 - exp(-n_N)))) / M;
end

% Fourier Machine theoretical BER calculation
function BER = BER_theoretical_Fourier(n_R, n_N, M)
    other_ports_clicks = 0;
    symbol = floor(M / 2);

    for m = 1:M

        if m ~= symbol
            I_mm_diff = n_R / (M * sin(pi * (symbol - m) / M) ^ 2);
            other_ports_clicks = other_ports_clicks + (1 - exp(- (I_mm_diff + n_N)));
        end

    end

    BER = (exp(- ((n_R * (2 * M ^ 2 + 1) / (3 * M)) + n_N)) + other_ports_clicks) / M;
end

