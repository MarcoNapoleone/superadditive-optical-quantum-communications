clear all
close all
nR=10;
M=[8 16 64]';

SNRdB = 0:0.1:20; % SNR in dB array of values [0, 1, 2, ..., 20]
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale
n_H=nR./SNR; % rumore Hadamard

for m=1:3
AUTO(m)=nR*(2*M(m).^2+1)./(3*M(m));
n_F(m,:)=nR./(SNR.*M(m)); % rumore Fourier

   for k=1:M(m)-1
CROSS(m,k)=nR/(M(m).*(sin(pi*k./M(m))).^2);
   end

% Hadamard (Green Machine) theoretical BER calculation
BER_H(m,:)=exp(- (M(m) * nR + n_H));

for k=1:M(m)-1
     BER_H(m,:)=BER_H(m,:)+(1 - exp(-n_H));
end
 BER_H(m,:)=BER_H(m,:)./M(m);

 % Fourier Machine theoretical BER calculation
BER_F(m,:)=exp(- (AUTO(m) + n_F(m,:)));
for k=1:M(m)-1
     BER_F(m,:)=BER_F(m,:)+(1 - exp(-CROSS(m,k)-n_F(m,:)));
end
 BER_F(m,:)=BER_F(m,:)./M(m);
 BER_PSK(m,:)=exp(-nR+n_H)/2+(1 - exp(-n_H))/2;
 SER_PSK(m,:)=1-(1-BER_PSK(m,:)).^log2(M(m));
end
figure (1)
semilogy(SNRdB,BER_F)
hold on
semilogy(SNRdB,BER_H, ':')
semilogy(SNRdB,SER_PSK,'--')

grid on
xlabel('SNR [dB]')
ylabel('SER')
legend ('Hadamard M=8','Hadamard M=16','Hadamard M=64', 'Fourier M=8','Fourier M=16','Fourier M=64','BPSK M=8','BPSK M=16','BPSK M=64')

