%% Simulation of the Green Machine and Fourier Machine
%
% This script simulates and compares the Bit Error Rate (BER) performance of 
% two different quantum communication machines: the Green Machine (Hadamard-based) 
% and the Fourier Machine. The simulation is carried out for different SNR values, 
% and both theoretical and Monte Carlo simulation results are presented.
%
% Machines:
%   - Green Machine: Based on Hadamard coding.
%   - Fourier Machine: Based on Fourier coding.
%
% Parameters:
%   - Planck constant, light speed, laser wavelength, signal power, baud rate.
%   - Monte Carlo simulations are used to estimate the BER under noisy conditions.
%
% The script performs the following steps:
%   1. Set up the physical parameters for the system.
%   2. Define the theoretical BER formulas for both machines.
%   3. Simulate the transmission over a noisy channel using Monte Carlo methods.
%   4. Plot the theoretical and simulated BER results for comparison.
%

clear;
close all;


%% Parameters
h = 6.62607015e-34;   % Planck constant [J·Hz^(-1)]
c = 299792458;        % Light speed [m/s]
lambda = 1550e-9;     % Wavelength of the laser [m]
nu = c / lambda;      % Frequency of the laser [Hz]
B = 1e8;              % Baud rate in Hz
P = 2e-13;            % Signal power in W
M = 32;               % Number of codewords Hadamard/Fourier
N_symbols = 8192;      % Number of transmission frames
m = floor(M / 2);     % Selected port for simulation


%% Channel parameters
SNRdB = 0:1:15; % SNR in dB array of values [0, 1, 2, ..., 20]
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale


%% Monte Carlo settings
N_simulations = 10e5;          % Number of Monte Carlo simulations
epsilon = 1e-4;       % Early stopping criterion for Monte Carlo simulation


%% Preallocating the symbols
% Randomly generate symbols to be transmitted from M possible codewords
rng(42);  % Set seed for reproducibility
frame = randi([0 M-1], 1, N_symbols); % message of length N_frames with values between 0 and M-1


%% Preallocate arrays for BER
BER_teo_H = zeros(1, length(SNR)); % Theoretical BER for Green Machine (Hadamard)
BER_teo_F = zeros(1, length(SNR)); % Theoretical BER for Fourier Machine
BER_sim_H = zeros(1, length(SNR)); % Simulated BER for Green Machine
BER_sim_F = zeros(1, length(SNR)); % Simulated BER for Fourier Machine


%% Photon number calculations for each machine
n_R_H = P / (h * nu * B) * log2(M); % Mean photon number for Green Machine
n_R_F = P / (h * nu * B) * log2(M); % Mean photon number for Fourier Machine

fprintf("The mean photon number n_R is: %f\n", n_R_H)


%% Theoretical BER Calculation Functions
% Hadamard (Green Machine) theoretical BER calculation
function BER = BER_theoretical_Hadamard(n_R, n_N, M)
    BER = (exp(- (M * n_R + n_N)) + ((M - 1) * (1 - exp(-n_N)))) / M;
end

% Fourier Machine theoretical BER calculation
function BER = BER_theoretical_Fourier(n_R, n_N, M)
    other_ports_clicks = 0;
    symbol = floor(M / 2);
    for m = 0:M-1
        if m ~= symbol
            I_mm_diff = n_R / (M * sin((pi * (symbol - m)) / M)^ 2);
            other_ports_clicks = other_ports_clicks + (1 - exp(- (I_mm_diff + n_N)));
        end
    end
    BER = (exp(- ((n_R * (2 * (M ^ 2) + 1) / (3 * M)) + n_N)) + other_ports_clicks) / M;
end


%% Simulation Functions
function errorCount = simulateHadamard(M, n_R, n_N, message, port)
    % Simulates the Hadamard Machine for one transmission frame
    errorCount = 0;
    for symbol = message
        if (symbol == port)  % Correct port
            if rand() < exp(- (M * n_R + n_N))        
                errorCount = errorCount + 1;
            end
                else  % Incorrect port
            if rand() < (1 - exp(-n_N))
                errorCount = errorCount + 1;
            end
        end
    end
end

function errorCount = simulateFourier(M, n_R, n_N, message, port)
    % Simulates the Fourier Machine for one transmission frame
    errorCount = 0;
    for symbol = message
        if (symbol == port)  % Correct port
            if rand() < exp(- ((n_R * (2 * (M ^ 2) + 1) / (3 * M)) + n_N))
                errorCount = errorCount + 1;
            end
        else  % Incorrect port
            I_mm_diff = n_R / (M * sin((pi * (port - symbol)) / M)^ 2);
            if rand() < (1 - exp(- (I_mm_diff + n_N)))
                errorCount = errorCount + 1;
            end
        end
    end
end


%% Loop over the SNR values
for idx = 1:length(SNR)
    snr = SNR(idx);

    % Calculate noise photons for Green and Fourier machines
    n_N_H = n_R_H / snr;         % Noise for Hadamard
    n_N_F = n_R_F / (snr * M);   % Noise for Fourier

    %% Theoretical BER (Green Machine and Fourier Machine)
    BER_teo_H(idx) = BER_theoretical_Hadamard(n_R_H, n_N_H, M);
    BER_teo_F(idx) = BER_theoretical_Fourier(n_R_F, n_N_F, M);

    %% Monte Carlo Simulation for Green Machine and Fourier Machine
    % Green Machine (Hadamard)
    simFunc_H = @() simulateHadamard(M, n_R_H, n_N_H, frame, m);

    % Fourier Machine
    simFunc_F = @() simulateFourier(M, n_R_F, n_N_F, frame, m);

    % Monte Carlo Simulator instances
    mc_H = MonteCarloSimulator(N_simulations, simFunc_H, epsilon);
    mc_F = MonteCarloSimulator(N_simulations, simFunc_F, epsilon);

    % Run the Monte Carlo simulations
    mc_H = mc_H.runSimulations();
    mc_F.verbose = true;
    mc_F = mc_F.runSimulations();

    % Store BER from Monte Carlo simulations
    BER_sim_H(idx) = mc_H.calculateMean() / N_symbols;
    BER_sim_F(idx) = mc_F.calculateMean() / N_symbols;
end




%% Plot Results
figure;
semilogy(SNRdB, BER_teo_H, 'b-.', 'LineWidth', 1);
hold on;
scatter(SNRdB, BER_sim_H , 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 1);
semilogy(SNRdB, BER_teo_F, 'r-.', 'LineWidth', 1);
scatter(SNRdB, BER_sim_F , 'rd', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 1);
hold off;

% Graph settings
xlabel('SNR (dB)');
ylabel('BER');
legend('Theoretical BER Hadamard', 'Simulated BER Hadamard', ...
       'Theoretical BER Fourier', 'Simulated BER Fourier');
title('Comparison of Theoretical vs Simulated BER for Hadamard and Fourier');
grid on;