clear;
close all;

%% Parameters
M = 256; % Fissato a 32
n_R_array = [ 0.005, 0.01, 0.05, 0.1, 1]; % Variazione di n_R

%% Channel parameters
SNRdB = 0:1:20; % SNR in dB array
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale

%% Preallocate arrays for BER
BER_teo_H = zeros(length(n_R_array), length(SNR)); % Theoretical BER for Hadamard (HM)
BER_teo_F = zeros(length(n_R_array), length(SNR)); % Theoretical BER for Fourier (FM)

%% Loop over n_R and SNR values
for nr_idx = 1:length(n_R_array)
    n_R = n_R_array(nr_idx); % Current n_R value
    fprintf("The mean photon number n_R is: %f\n", n_R)
    
    for idx = 1:length(SNR)
        snr = SNR(idx);

        % Calculate noise photons for Green and Fourier machines
        n_N_H = n_R / snr; % Noise for Hadamard
        n_N_F = n_R / (snr * M); % Noise for Fourier

        %% Theoretical BER (Hadamard and Fourier Machines)
        BER_teo_H(nr_idx, idx) = BER_theoretical_Hadamard(n_R, n_N_H, M);
        BER_teo_F(nr_idx, idx) = BER_theoretical_Fourier(n_R, n_N_F, M);
    end
end

%% Plot the ratio for each n_R
figure;
styles = {'-', '--', '-.', ':', '-'}; % Line styles for different n_R values
colors = lines(length(n_R_array)); % Generate distinct colors

for nr_idx = 1:length(n_R_array)
    ratio_HF = BER_teo_H(nr_idx, :) ./ BER_teo_F(nr_idx, :); % Ratio between Hadamard and Fourier BER
    plot(SNRdB, ratio_HF, 'LineStyle', styles{nr_idx}, 'Color', colors(nr_idx, :), 'LineWidth', 1.5, ...
        'DisplayName', ['n_R = ' num2str(n_R_array(nr_idx))]);
    hold on;
end

%dash hor line at 1
grayColor = [.5 .5 .5];
line([min(SNRdB), max(SNRdB)], [1, 1], 'Color', grayColor, 'LineStyle', ':', 'LineWidth', 1, DisplayName='Ratio = 1');
xlabel('SNR (dB)');
ylabel('Ratio (Hadamard/Fourier)');
title('Ratio of SER for (Hadamard/Fourier) Machines at varying n_R');
grid on;
legend();

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
