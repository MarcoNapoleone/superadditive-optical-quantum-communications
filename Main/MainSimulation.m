addpath('../MonteCarlo')

clear;
close all;

%% Parameters
h = 6.62607015e-34; % Planck constant [J·Hz^(-1)]
c = 299792458; % Light speed [m/s]
lambda = 1550e-9; % Wavelength of the laser [m]
nu = c / lambda; % Frequency of the laser [Hz]
B = 10e8; % Baud rate in Hz
P = 1e-12; % Signal power in W
M = 128; % Number of codewords Hadamard/Fourier
N_symbols = 8192; % Number of transmission frames
m = floor(M / 2); % Selected port for simulation

%% Channel parameters
SNRdB = 0:1:20; % SNR in dB array of values [0, 1, 2, ..., 20]
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale

%% Monte Carlo settings
N_simulations = 10e7; % Number of Monte Carlo simulations
epsilon = 1e-5; % Early stopping criterion for Monte Carlo simulation

%% Preallocating the symbols
% Randomly generate symbols to be transmitted from M possible codewords
rng(42); % Set seed for reproducibility
frame = randi([0 M - 1], 1, N_symbols); % message of length N_frames with values between 0 and M-1

%% Preallocate arrays for BER and standard deviation (std)
BER_teo_H = zeros(1, length(SNR)); % Theoretical BER for Green Machine (Hadamard)
BER_teo_F = zeros(1, length(SNR)); % Theoretical BER for Fourier Machine

BER_sim_H = zeros(1, length(SNR)); % Simulated BER for Green Machine
BER_sim_F = zeros(1, length(SNR)); % Simulated BER for Fourier Machine

BER_std_H = zeros(1, length(SNR)); % Standard deviation for Green Machine
BER_std_F = zeros(1, length(SNR)); % Standard deviation for Fourier Machine

BER_Q_H = zeros(1, length(SNR)); % Quartile for Green Machine
BER_Q_F = zeros(1, length(SNR)); % Quartile for Fourier Machine

%% Photon number calculations for each machine
n_R_H = P / (h * nu * B) * log2(M); % Mean photon number for Green Machine
n_R_F = P / (h * nu * B) * log2(M); % Mean photon number for Fourier Machine

fprintf("The mean photon number n_R is: %f\n", n_R_H)

%% Loop over the SNR values
for idx = 1:length(SNR)
    snr = SNR(idx);

    % Calculate noise photons for Green and Fourier machines
    n_N_H = n_R_H / snr; % Noise for Hadamard
    n_N_F = n_R_F / (snr * M); % Noise for Fourier

    %% Theoretical BER (Green Machine and Fourier Machine)
    BER_teo_H(idx) = BER_theoretical_Hadamard(n_R_H, n_N_H, M);
    BER_teo_F(idx) = BER_theoretical_Fourier(n_R_F, n_N_F, M);

    %% Monte Carlo Simulation for Green Machine and Fourier Machine
    % Green Machine (Hadamard)
    simFunc_H = @() simulateHadamard(n_R_H, n_N_H, M, frame, m);

    % Fourier Machine
    simFunc_F = @() simulateFourier(n_R_F, n_N_F, M, frame, m);

    % Monte Carlo Simulator instances
    mc_H = MonteCarloSimulator(N_simulations, simFunc_H, epsilon);
    mc_F = MonteCarloSimulator(N_simulations, simFunc_F, epsilon);

    % Run the Monte Carlo simulations and gather statistics
    mc_H = mc_H.runSimulations();
    mc_F.verbose = true;
    mc_F = mc_F.runSimulations();

    % Store BER and standard deviation from Monte Carlo simulations
    [mean_H, std_H, q_H] = mc_H.calculateResult();
    [mean_F, std_F, q_F] = mc_F.calculateResult();

    BER_sim_H(idx) = mean_H / N_symbols;
    BER_sim_F(idx) = mean_F / N_symbols;

end

%% Plot Results with Error Bars
figure;

colors = {'#098bf8', '#fe5f55'}; % Colors for Green and Fourier Machines

% Theoretical curves
semilogy(SNRdB, BER_teo_H, '-.', 'LineWidth', 1.5, 'Color', colors{1}, 'DisplayName', "Theoretical Hadamard");
hold on;
scatter(SNRdB, BER_sim_H , 'o', 'MarkerFaceColor', colors{1}, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', colors{1}, 'MarkerEdgeAlpha', 1, 'DisplayName', "Numerical Hadamard");

semilogy(SNRdB, BER_teo_F, '-', 'LineWidth', 1.5, 'Color', colors{2}, 'DisplayName', "Theoretical Fourier");
scatter(SNRdB, BER_sim_F , 'd', 'MarkerFaceColor', colors{2}, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', colors{2}, 'MarkerEdgeAlpha', 1, 'DisplayName', "Numerical Fourier");
hold off;

% Graph settings
legend();
xlabel('SNR (dB)');
ylabel('SER');
grid on;


%export pdf
exportgraphics(gcf, 'output/BERvsSNR_MonteCarlo.pdf')

%% Theoretical BER Calculation Functions
% Hadamard (Green Machine) theoretical BER calculation
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


%% Simulation Functions
function errorCount = simulateHadamard(n_R, n_N, M, message, port)
    % Simulates the Hadamard Machine for one transmission frame
    errorCount = 0;

    for symbol = message

        if (symbol == port) % Correct port

            if rand(1) < exp(- (M * n_R + n_N))
                errorCount = errorCount + 1;
            end

        else % Incorrect port

            if rand(1) < (1 - exp(-n_N))
                errorCount = errorCount + 1;
            end

        end

    end

end

function errorCount = simulateFourier(n_R, n_N, M, message, port)
    % Simulates the Fourier Machine for one transmission frame
    errorCount = 0;

    for symbol = message

        if (symbol == port) % Correct port

            if rand(1) < exp(- ((n_R * (2 * M ^ 2 + 1) / (3 * M)) + n_N))
                errorCount = errorCount + 1;
            end

        else % Incorrect port
            I_mm_diff = n_R / (M * sin((pi * (port - symbol)) / M) ^ 2);

            if rand(1) < (1 - exp(- (I_mm_diff + n_N)))
                errorCount = errorCount + 1;
            end

        end

    end

end
