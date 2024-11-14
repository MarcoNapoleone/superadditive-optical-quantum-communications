# Superadditive Optical Quantum Communications

## Overview
This project explores the theoretical foundations and challenges of integrating quantum superadditive technologies into existing optical networks, evaluating their performance in a noisy channel. The implementation is based on the concepts discussed in the paper "Perspectives on Superadditive Optical Quantum Communications" by Gabriella Cincotti, Marco Napoleone, and Matteo Rosati.

## Key Concepts
1. **Superadditive Quantum Communications**: Utilizes joint detection receivers (JDR) to enhance system capacity beyond classical limits.
2. **Hadamard machine (Green machine)**: Structured optical receivers using Hadamard codewords to improve spectral efficiency.
3. **Fourier machine**: Novel JDR based on Fourier codes, offering advantages in low-power optical settings.

## Implementation Details
- **Languages**: MATLAB
- **Frameworks**: None

## Project Structure
- `Main/`: Contains the source code for the implementation.
- `MonteCarlo/`: Monte Carlo class for simulating the performance of the system.
- `Hadamard/`: Implementation of the Hadamard codes.
- `OOK/`: Implementation of the on-off keying modulation.
- `README.md`: Project documentation.

## Results
The project evaluates the performance of the superadditive quantum communication system in a noisy channel, comparing the Hadamard and Fourier machines.
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Images in a Row</title>
    <style>
        .image-row img {
            width: 200px;
            display: inline-block;
        }
    </style>
</head>
<body>
    <div class="image-row">
        <img src="Main/output/BERvsSNR_M-1.png" alt="BER vs SNR for M=1">
        <img src="Main/output/BERvsSNR_nR_y-1.png" alt="BER vs SNR for nR_y=1">
        <img src="Main/output/BERvsSNR_MonteCarlo-1.png" alt="BER vs SNR for M=1 (Monte Carlo)">
    </div>
</body>
</html>

## Contact
For any questions or issues, please contact Gabriella Cincotti at gabriella.cincotti@uniroma3.it.