clear;
close all;

%% Parameters
M = 128; % Fixed M value for Hadamard/Fourier Machine

% Array of mean photon numbers
n_R = [0.01,0.1, 1]

%% Channel parameters
SNRdB = 0:1:20; % SNR in dB array
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale

%% Preallocate arrays for BER
BER_teo_H = zeros(length(n_R), length(SNR)); % Theoretical BER for Hadamard (HM)
BER_teo_F = zeros(length(n_R), length(SNR)); % Theoretical BER for Fourier (FM)

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
    end

end

%% Plot Results for each n_R
figure;
styles = {'-', '-.', '--', ':'}; % Line styles for different M
colors = {'#098bf8', '#fe5f55'}; % Colors for Green and Fourier Machines

for nR_idx = 1:length(n_R)
    current_n_R = n_R(nR_idx);

    % Plot Hadamard results for current n_R
    plot(SNRdB, BER_teo_H(nR_idx, :), 'Color', colors{1}, 'LineStyle', styles{nR_idx}, 'LineWidth', 1.5, displayname = 'Hadamard n_R=' + string(current_n_R));
    hold on;

    % Plot Fourier results for current n_R
    semilogy(SNRdB, BER_teo_F(nR_idx, :), 'Color', colors{2}, 'LineStyle', styles{nR_idx}, 'LineWidth', 1.5, displayname = 'Fourier n_R=' + string(current_n_R));
end

% Graph settings
xlabel('SNR (dB)');
ylabel('SER');


legend();
grid on;

%export pdf
exportgraphics(gcf, 'output/BERvsSNR_nR.pdf')



%% Theoretical BER Calculation Functions
% Hadamard (Green Machine) theoretical BER calculation
function BER = BER_theoretical_Hadamard(n_R, n_N, M)
    BER = (exp(- (M * n_R + n_N)) + ((M - 1) * (1 - exp(-n_N)))) / M;
end

% Fourier Machine theoretical BER calculation
function BER = BER_theoretical_Fourier(n_R, n_N, M)
    other_ports_clicks = 0;
    symbol = floor(M / 2);

    for m = 0:M - 1

        if m ~= symbol
            I_mm_diff = n_R / (M * sin(pi * (symbol - m) / M) ^ 2);
            other_ports_clicks = other_ports_clicks + (1 - exp(- (I_mm_diff + n_N)));
        end

    end

    BER = (exp(- ((n_R * (2 * M ^ 2 + 1) / (3 * M)) + n_N)) + other_ports_clicks) / M;
end
