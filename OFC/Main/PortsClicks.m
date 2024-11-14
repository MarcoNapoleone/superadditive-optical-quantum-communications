%% Visualizzazione delle probabilità di click per diverse SNR
%
% Questo script visualizza la probabilità di click per ciascuna porta delle
% due macchine (Green Machine e Fourier Machine) su quattro subplot, ciascuno con un diverso valore di SNR.
%

clear;
close all;

%% Parametri
M = 32;               % Numero di porte (codewords Hadamard/Fourier)
m = floor(M / 2);     % La porta selezionata (porta corretta)
SNRdB_list = [5, 10, 15, 20]; % Diversi valori di SNR in dB

%% Numero medio di fotoni
n_R = 0.1;

%% Preparazione della figura
figure;

%% Loop sui diversi valori di SNR
for i = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i);          % SNR in dB
    SNR = 10^(SNRdB / 10);          % Conversione in scala lineare
    n_N_H = n_R / SNR;            % Fotoni di rumore per Green Machine
    n_N_F = n_R / (SNR * M);      % Fotoni di rumore per Fourier Machine

    %% Probabilità di click per ciascuna porta
    clickProb_H = zeros(1, M); % Probabilità per Green Machine (Hadamard)
    clickProb_F = zeros(1, M); % Probabilità per Fourier Machine

    %% Calcolo della probabilità di click per la Green Machine (Hadamard)
    for port = 1:M
        if port == m
            % Porta corretta
            clickProb_H(port) = exp(- (M * n_R + M*n_N_H));
        else
            % Porte sbagliate
            clickProb_H(port) = 1 - exp(-M*n_N_H);
        end
    end

    %% Calcolo della probabilità di click per la Fourier Machine
    for port = 1:M
        if port == m
            % Porta corretta
            clickProb_F(port) = exp(- (n_R * (2 * (M^2) + 1) / (3 * M)) - M*n_N_F);
        else
            % Porte sbagliate
            I_mm_diff = n_R / (M * sin((pi * (port - m)) / M)^2);
            clickProb_F(port) = 1 - exp(- (I_mm_diff + M*n_N_F));
        end
    end



    %% Plot per la Green Machine
    subplot(4, 2, 2 * i - 1); % Subplot dispari per Green Machine
    hold on;
    % Disegna la barra pastello (complementare)
    bar(1:M, ones(1, M), 'FaceColor', [1, 0.75, 0.75], 'EdgeColor', 'none'); % Sfondo pastello fino a 1
    % Disegna la barra rossa in foreground con il valore reale
    bar(1:M, clickProb_H, 'FaceColor', [0.5, 0, 0], 'EdgeColor', 'none'); % Probabilità reale in rosso
    title(['Green Machine (SNR = ', num2str(SNRdB), ' dB)']);
    xlabel('Porta');
    ylabel('Probabilità di Click');
    xticks(1:M);
    grid on;

    %% Plot per la Fourier Machine
    subplot(4, 2, 2 * i); % Subplot pari per Fourier Machine
    hold on;
    % Disegna la barra pastello (complementare) blu 
    bar(1:M, ones(1, M), 'FaceColor', [0.75, 0.75, 1], 'EdgeColor', 'none'); % Sfondo pastello fino a 1
    % Disegna la barra rossa in foreground con il valore reale
    bar(1:M, clickProb_F, 'FaceColor', [0, 0, 0.5], 'EdgeColor', 'none'); % Probabilità reale in blu
    title(['Fourier Machine (SNR = ', num2str(SNRdB), ' dB)']);
    xlabel('Porta');
    ylabel('Probabilità di Click');
    xticks(1:M);
    grid on;

end

%% Impostazioni grafiche generali
sgtitle('Probabilità di Click per Green Machine e Fourier Machine a Diversi Valori di SNR');

