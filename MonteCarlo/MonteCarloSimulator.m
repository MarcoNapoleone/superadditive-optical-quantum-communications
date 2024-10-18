%% MonteCarloSimulator Class
% 
% This class performs Monte Carlo simulations with an option for early stopping
% based on a predefined threshold (epsilon). The user can define a simulation 
% function, set the number of simulations, and analyze the results using built-in
% methods for calculating the mean and standard deviation.
%
% Properties:
%   - numSimulations:  The maximum number of Monte Carlo simulations to run.
%   - simulationFunc:  A function handle that defines the simulation logic.
%   - epsilon:         A threshold for early stopping, based on the stability of 
%                      the running mean of the simulation results.
%   - results:         An array storing the results of each simulation run.
%   - verbose:         A flag to indicate whether early stopping messages should be printed.
%
% Methods:
%   - MonteCarloSimulator:  Class constructor to initialize properties.
%   - runSimulations:       Runs the Monte Carlo simulations with early stopping.
%   - calculateMean:        Calculates the mean of the simulation results.
%   - calculateStd:         Calculates the standard deviation of the simulation results.
%   - plotResults:          Plots the histogram of the simulation results.

classdef MonteCarloSimulator
    properties
        % numSimulations: Maximum number of simulations to run.
        numSimulations
        
        % simulationFunc: Handle to the simulation function, provided by the user.
        simulationFunc
        
        % epsilon: Threshold for early stopping when the mean stabilizes.
        epsilon
        
        % results: Array to store the results of each simulation run.
        results
        
        % verbose: Flag to control the printing of early stopping messages.
        verbose = false
    end
    
    methods
        % Constructor for the MonteCarloSimulator class
        % 
        % Parameters:
        %   numSimulations (int): Number of simulations to perform.
        %   simulationFunc (function_handle): Function handle for the simulation logic.
        %   epsilon (double): Early stopping threshold based on the mean difference.
        %
        % Returns:
        %   obj: An instance of the MonteCarloSimulator class.
        function obj = MonteCarloSimulator(numSimulations, simulationFunc, epsilon)
            if nargin == 3  % Check if all arguments are provided
                obj.numSimulations = numSimulations;
                obj.simulationFunc = simulationFunc;
                obj.epsilon = epsilon;
                obj.results = zeros(1, numSimulations);  % Pre-allocate results array
            else
                error('MonteCarloSimulator requires 3 inputs: numSimulations, simulationFunc, and epsilon.');
            end
        end
        
        % Runs the Monte Carlo simulations with early stopping.
        %
        % The simulation stops early if the mean of the results converges,
        % i.e., the change in mean is smaller than the epsilon threshold.
        %
        % Parameters:
        %   verbose (bool, optional): Set to true to print early stopping information. Default is false.
        function obj = runSimulations(obj)
            previousMean = 0;  % Initialize the previous mean
            
            for i = 1:obj.numSimulations
                obj.results(i) = obj.simulationFunc();  % Run the simulation
                
                % Calculate the current mean up to this point
                currentMean = mean(obj.results(1:i));
                
                % Check for early stopping if the difference in means is below epsilon
                if i > 1 && abs(currentMean - previousMean) < obj.epsilon
                    if obj.verbose
                        fprintf('Early stopping criterion met at iteration %d, ', i);
                        fprintf('Mean calculated: %.4f\n', currentMean);
                    end
                    obj.results = obj.results(1:i);  % Trim results to the number of iterations run
                    break;
                end
                
                % Update the previous mean
                previousMean = currentMean;
            end
        end
        
        % Calculates the mean of the simulation results.
        
        function [meanResult, stdResult, quartiles] = calculateResult(obj)
            meanResult = mean(obj.results);
            stdResult = std(obj.results);
            quartiles = quantile(obj.results, [0.25, 0.75]);
        end

        % Plots the results of the simulations as a histogram.
        %
        % This method creates a probability histogram of the simulation results.
        function plotResults(obj)
            figure;
            histogram(obj.results, 'Normalization', 'probability');
            xlabel('Simulation Results');
            ylabel('Probability');
            title('Monte Carlo Simulation Results');
        end
    end
end
