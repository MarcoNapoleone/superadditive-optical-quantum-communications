clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols

%% Tx parameters
modfor = 'ook';         % modulation format    (OnOff keying)

%% Channel parameters
SNRdB = 5:1:14;         % signal to noise ratio [dB]        

%% Rx parameters
rx.type = 'bin';        % get binary patterns

%% Monte Carlo BER parameters

mc.stop = [0.1 68];     % stop at relative error 0.1 with confidence 68%
mc.maxsamp = 1e8;       % maximum number of samples for Monte Carlo

%% Init
inigstate(Nsymb);       % initialize global variable

avgber = zeros(1,length(SNRdB));    % average bit-error rate
stdber = zeros(1,length(SNRdB));    % std of avgber (accuracy of estimation)

%% Discrete-time AWGN channel

rng(1);
[pat, patbin] = pattern2(Nsymb,'rand',struct('format',modfor)); % Generate random binary pattern for OOK
ak = pat2samp(pat,modfor); % Convert pattern to OOK samples (0 or 2 for OOK) avg energy = 1

for k=1:length(SNRdB)

    %% complex noise variance  (SNRdb = 10*log10(1/sigma2), assuming avg energy of ak = 1)
    sigma2 = 1/10^(SNRdB(k)/10); 
    cond = true; 
    while cond % when cond=false the desired accuracy for avgber has been reached
        
        % Generate AWGN noise
        nk = sqrt(sigma2)*(randn(Nsymb,1)+1i*randn(Nsymb,1));        
        rk = ak + nk; % AWGN
        %fprintf('rk = %d, ak = %d, nk = %d\n', rk, ak, nk);

        % Decision rule for OOK: 
        patbinhat = double(real(rk) >= 1); % Binary decisions at Rx (OOK detection)
        
        % Count bit-errors, i.e., estimate the BER
        [cond,avgber(k),nruns,stdber(k)] = meanestimate(patbin ~= patbinhat,mc);
        
        % Display info
        fprintf('   SNR = %4.1f [dB]. Bits observed = %6d, <BER> = %.2e, std(BER)/<BER> = %.3f.\n',...
            SNRdB(k),nruns,avgber(k),stdber(k)/avgber(k));
    end
    fprintf('\n')
end

%% Plot

sigma2 = 1./(10.^(SNRdB/10));  % complex noise variance
bertheory = 0.5*erfc(1./sqrt(sigma2*2)); % theoretical BER for OOK for [0, 2] levels
semilogy(SNRdB,avgber,'o',SNRdB,bertheory,'k--')
grid on
xlabel('SNR [dB]')
ylabel('BER')
legend('estimation','theory')
