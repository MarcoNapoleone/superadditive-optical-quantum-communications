clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]. (bitrate)
modfor = 'ook';         % modulation format    (OnOff keying)

%% Channel parameters
SNRdB = [3, 5, 7, 9];  % different SNR values [dB]        

%% Rx parameters
rx.type = 'bin';        % get binary patterns

%% Monte Carlo BER parameters

mc.stop = [0.1 68];     % stop at relative error 0.1 with confidence 68%
mc.maxsamp = 1e8;       % maximum number of samples for Monte Carlo

%% Init
inigstate(Nsymb);       % initialize global variable

avgber = zeros(1,length(SNRdB));    % average bit-error rate
stdber = zeros(1,length(SNRdB));    % std of avgber (accuracy of estimation)

%% Create subplots
figure;

% Generate random binary pattern for OOK
[pat, patbin] = pattern2(Nsymb,'rand',struct('format',modfor));
ak = pat2samp(pat,modfor); % Convert pattern to OOK samples (0 or 2 for OOK) avg energy = 1


for i = 1:length(SNRdB)
    % Extract current SNR value
    currentSNRdB = SNRdB(i);
    sigma2 = 1/10^(currentSNRdB/10); % complex noise variance

    cond = true; 
    while cond % whe
        % Generate AWGN noise
        nk = sqrt(sigma2)*(randn(Nsymb,1)+1i*randn(Nsymb,1));        
        rk = ak + nk; % AWGN

        % Decision rule for OOK: 
        patbinhat = double(real(rk) >= 1); % Binary decisions at Rx (OOK detection)
        
        % Count bit-errors, i.e., estimate the BER
        [cond,avgber(i),nruns,stdber(i)] = meanestimate(patbin ~= patbinhat,mc);

        fprintf('   SNR = %4.1f [dB]. Bits observed = %6d, <BER> = %.2e, std(BER)/<BER> = %.3f.\n',...
            currentSNRdB,nruns,avgber(i),stdber(i)/avgber(i));
    end

    % Create subplot
    subplot(ceil(length(SNRdB)/2), 2, i);
    
    % Plot symbols
    scatter(real(rk(ak == 0)), imag(rk(ak == 0)), 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); % Red for symbol 0
    hold on;
    scatter(real(rk(ak == 2)), imag(rk(ak == 2)), 'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); % Blue for symbol 2
    
    % Add the decision threshold line (vertical line at Re = 1)
    xline(1, '--r', 'Threshold');

    % Plot formatting
    title(sprintf('SNR = %d dB\nBER = %.2e', currentSNRdB, avgber(i)));
    xlabel('Re(rk) (Real part)');
    ylabel('Im(rk) (Imaginary part)');
    grid on;
    legend('Symbol 0', 'Symbol 2', 'Threshold (Re = 1)');
    hold off;
end
