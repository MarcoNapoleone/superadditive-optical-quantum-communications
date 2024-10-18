clear all
close all
M=16;
nR=[0.001 0.1 1]';

SNRdB = 0:1:15; % SNR in dB array of values [0, 1, 2, ..., 20]
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale
for m=1:3
n_H=nR(m)./SNR; % rumore Hadamard

AUTO(m)=nR(m)*(2*M.^2+1)./(3*M);
n_F(m,:)=nR(m)./(SNR.*M); % rumore Fourier

   for k=1:M-1
CROSS(m,k)=nR(m)/(M.*(sin(pi*k./M)).^2);
   end

% Hadamard (Green Machine) theoretical BER calculation
BER_H(m,:)=exp(- (M * nR(m) + n_H));

for k=1:M-1
     BER_H(m,:)=BER_H(m,:)+(1 - exp(-n_H));
end
 BER_H(m,:)=BER_H(m,:)./M;

 % Fourier Machine theoretical BER calculation
BER_F(m,:)=exp(- (AUTO(m) + n_F(m,:)));
for k=1:M-1
     BER_F(m,:)=BER_F(m,:)+(1 - exp(-CROSS(m,k)-n_F(m,:)));
end
 BER_F(m,:)=BER_F(m,:)./M;
end
figure (1)
semilogy(SNRdB,BER_F)
legend('Fourier')
hold on
semilogy(SNRdB,BER_H, ':')
grid on
xlabel('SNR [dB]')
ylabel('BER')
legend('Hadamard')

