%% Monte Carlo Simulation Example: Rolling Two Dice
%
% This script demonstrates how to use the MonteCarloSimulator class to perform
% a Monte Carlo simulation of rolling two dice and summing the result. The simulation
% is run for a defined number of iterations, with an option for early stopping 
% if the mean of the results converges based on a specified threshold (epsilon).
%
% The key steps are:
%   1. Defining a simulation function that simulates rolling two dice.
%   2. Running the Monte Carlo simulation with early stopping.
%   3. Calculating and displaying the mean and standard deviation of the results.
%   4. Plotting a histogram of the simulation outcomes.


% Clear workspace and close all figures
clear;
close all;

% Step 1: Define the number of simulations to perform
numSimulations = 10e8;  % Large number for a robust Monte Carlo simulation

% Step 2: Define the simulation function
% The simulation function simulates the sum of rolling two 6-sided dice.
simulationFunc = @() randi([1, 6]) + randi([1, 6]);

% Step 3: Define the early stopping criterion (epsilon)
% The simulation will stop early if the change in mean is smaller than this value.
epsilon = 10e-6;

% Step 4: Create a Monte Carlo Simulator instance
% Initialize the MonteCarloSimulator with the number of simulations, simulation function, and epsilon.
mc = MonteCarloSimulator(numSimulations, simulationFunc, epsilon);

% Step 5: Run the simulations with early stopping
% The simulation will stop early if the mean stabilizes.
mc = mc.runSimulations();

% Step 6: Calculate the mean and standard deviation of the simulation results
[meanResult,stdResult, ~] = mc.calculateResult();  % Calculate the mean of the results

% Display the calculated mean and standard deviation
fprintf('Mean of the results: %.4f\n', meanResult);
fprintf('Standard deviation of the results: %.4f\n', stdResult);

% Step 7: Plot the histogram of the simulation results
% This will display the probability distribution of the outcomes.
mc.plotResults();
