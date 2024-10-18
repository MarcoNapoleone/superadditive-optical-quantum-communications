clear;
close all;

%% Parameters

%M is array from 8, 16, 32, 64, 128
p = 4:1:7; 
M = 2 .^ p; % Number of codewords Hadamard/Fourier

% Photon number calculations for each machine
n_R = 0.1


%% Channel parameters
SNRdB = 0:1:15; % SNR in dB array
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale


%% Preallocate arrays for BER
BER_teo_H = zeros(length(M), length(SNR)); % Theoretical BER for Hadamard (HM)
BER_teo_F = zeros(length(M), length(SNR)); % Theoretical BER for Fourier (FM)

fprintf("The mean photon number n_R is: %f\n", n_R)


%% Loop over each M and SNR values
for m_idx = 1:length(M)
    current_M = M(m_idx); % Select current M
    for idx = 1:length(SNR)
        snr = SNR(idx);

        % Calculate noise photons for Green and Fourier machines
        n_N_H = n_R / snr; % Noise for Hadamard
        n_N_F = n_R / (snr * current_M); % Noise for Fourier

        %% Theoretical BER (Hadamard and Fourier Machines)
        BER_teo_H(m_idx, idx) = BER_theoretical_Hadamard(n_R, n_N_H, current_M);
        BER_teo_F(m_idx, idx) = BER_theoretical_Fourier(n_R, n_N_F, current_M);
    end
end


%% Plot Results for each M
figure;
styles = {'-', '--', '-.', ':'}; % Line styles for different M
colors = {'#098bf8', '#fe5f55'}; % Colors for Green and Fourier Machines

% Blue color for Hadamard, Red color for Fourier

for m_idx = 1:length(M)
    current_M = M(m_idx);
        
    semilogy(SNRdB, BER_teo_H(m_idx, :), 'Color', colors{1}, 'LineStyle', styles{m_idx}, 'LineWidth', 2, 'DisplayName', 'Hadamard M=' + string(current_M));
    hold on;
    semilogy(SNRdB, BER_teo_F(m_idx, :), 'Color', colors{2}, 'LineStyle', styles{m_idx}, 'LineWidth', 2, 'DisplayName', 'Fourier M=' + string(current_M));
    
end

% Graph settings
xlabel('SNR (dB)');
ylabel('SER');
%leged on bottom left
legend('Location', 'southwest');
grid on;

%export pdf
exportgraphics(gcf, 'output/BERvsSNR_M.pdf')



%% Theoretical BER Calculation Functions
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

