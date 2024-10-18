

nR = 0.7802
M = 32;
AUTO = nR * (2 * M ^ 2 + 1) / (3 * M);

for k = 1:M - 1
    CROSS(k) = nR / (M * (sin(pi * k / M)) ^ 2);
end

SNRdB = 0:1:15; % SNR in dB array of values [0, 1, 2, ..., 20]
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale
n_H = nR ./ SNR; % rumore Hadamard
n_F = nR ./ (SNR * M); % rumore Fourier

% Hadamard (Green Machine) theoretical BER calculation
BER_H = exp(- (M * nR + n_H));

for k = 1:M - 1
    BER_H = BER_H + (1 - exp(-n_H));
end

BER_H = BER_H / M;

% Fourier Machine theoretical BER calculation
BER_F = exp(- (AUTO + n_F));

for k = 1:M - 1
    BER_F = BER_F + (1 - exp(-CROSS(k) - n_F));
end

BER_F = BER_F / M;
figure (1)
semilogy(SNRdB, BER_F, '-r', SNRdB, BER_H, '-k')
grid on
xlabel('SNR [dB]')
ylabel('BER')
legend('Fourier', 'Hadamard')
