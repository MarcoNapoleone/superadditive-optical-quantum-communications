# Simulation of the Green Machine and Fourier Machine

This document explains the simulation and comparison of the Bit Error Rate (BER) performance of two quantum communication systems: the **Green Machine** (based on Hadamard coding) and the **Fourier Machine** (based on Fourier coding). Both theoretical and Monte Carlo simulation results are presented for various Signal-to-Noise Ratio (SNR) values.

---

## 1. Introduction

Quantum communication systems require efficient encoding schemes for transmitting information over noisy channels. The **Green Machine** uses Hadamard coding, and the **Fourier Machine** uses Fourier coding. This simulation aims to analyze and compare their BER performance, which is critical for determining system reliability.

### Machines:
- **Green Machine**: Implements Hadamard coding for quantum transmission.
- **Fourier Machine**: Implements Fourier coding, which distributes the signal across frequency modes.

### Parameters:
- **Planck constant (h)**: Fundamental constant used in photon energy calculation.
- **Light speed (c)**: The speed of light in vacuum.
- **Wavelength (λ)**: Wavelength of the laser used in the system.
- **Signal power (P)**: Transmitted signal power.
- **Baud rate (B)**: Transmission rate.
- **Monte Carlo simulations**: Used to estimate BER performance under noisy conditions.

---

## 2. Photon Number Calculation

The number of photons per symbol, $n_R$, is a key factor determining the BER in a quantum communication system. It is calculated as:

$$
n_R = \frac{P}{h \cdot \nu \cdot B} \cdot \log_2(M)
$$

Where:
- $P$ is the signal power.
- $h$ is Planck's constant.
- $\nu = \frac{c}{\lambda}$ is the frequency of the laser.
- $B$ is the baud rate.
- $M$ is the number of codewords used (either Hadamard or Fourier).

The **mean photon number** is identical for both the Green and Fourier machines.

---

## 3. Theoretical BER Calculation

### 3.1 Hadamard Coding (Green Machine)
The theoretical BER for Hadamard coding can be calculated as:

$$
\text{BER}_{\text{Hadamard}} = \frac{\exp\left( -(M \cdot n_R + n\_N) \right) + (M - 1) \cdot \left(1 - \exp\left(-n\_N\right)\right)}{M}
$$

Where:
- $M$ is the number of codewords.
- $n_R$ is the mean photon number.
- $n_N$ is the number of noise photons.

### 3.2 Fourier Coding (Fourier Machine)
The theoretical BER for Fourier coding is given by:

$$
\text{BER}_{\text{Fourier}} = \frac{\exp \left( -\left( \frac{n_R \cdot (2 \cdot M^2 + 1)}{3 \cdot M} \right) + n_N \right) + \sum_{m=1}^{M-1} \left( 1 - \exp \left( - \left( I^{mm'} + n_N \right) \right) \right)}{M}
$$

Where $I_{mm\_diff}$ is the interference between different ports and is calculated as:

$$
I^{mm'} = \frac{n_R}{M \cdot \sin^2 \left( \frac{\pi \cdot (m - m'
)}{M} \right)}
$$

Here, $n_N$ is the noise photons for the Fourier Machine.

---

## 4. Monte Carlo Simulation

### 4.1 Green Machine (Hadamard) Simulation
For each symbol transmission, the simulation compares the transmitted symbol (codeword) to the selected port. If the correct port is selected, the probability of error is determined by:

$$
P_{\text{error}} = \exp \left( -(M \cdot n_R + n_N) \right)
$$

If the wrong port is selected, the probability of error is given by:

$$
P_{\text{error}} = 1 - \exp(-n_N)
$$

### 4.2 Fourier Machine Simulation
In the Fourier Machine, the error calculation depends on whether the correct port is selected or not. If the correct port is selected, the error probability is:

$$
P_{\text{error}} = \exp \left( -\left( \frac{n_R \cdot (2 \cdot M^2 + 1)}{3 \cdot M} + n_N \right) \right)
$$

For the incorrect port, the error probability is based on the interference between ports:

$$
P_{\text{error}} = 1 - \exp \left( - \left( I^{mm'} + n_N \right) \right)
$$

---

## 
