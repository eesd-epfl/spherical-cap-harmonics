%*************************************************************************%
% All rights reserved (C) to the authors: Mahmoud SHAQFA, Gary CHOI       %
%		 and Katrin BEYER                                                 %
%                                                                         %
% M. Shaqfa Contact:                                                      %
% Earthquake Engineering and Structural Dynamics Laboratory (EESD),       %
% School of Architecture, Civil and Environmental Engineering (ENAC),     %
% Ecole polytechnique federale de Lausanne (EPFL),                        %
% CH-1015 Lausanne, Switzerland.                                          %
%               Tel.: +41 21 69 33297                                     %
%               Email: mahmoud.shaqfa@epfl.ch                             %
%                                                                         %
% G. Choi Contact:                                                        %
% Department of Mathematics, Massachusetts Institute of Technology (MIT)  %
% Cambridge, MA, USA                                                      %
%               Email: ptchoi@mit.edu                                     %
%                                                                         %
% K. Beyer Contact:                                                       %
%               Email: katrin.beyer@epfl.ch                               %
%*************************************************************************%
% This code includes implementations for:                                 %
%				- Spherical Cap Harmonics (SCH)                           %
%				- Spherical Harmonics (SH)                                %
%				- HemiSpherical Harmonics (HSH)                           %
% This code is part of the paper: "Spherical Cap Harmonic Analysis(SCHA)..%
%	 for Characterising the Morphology of Rough Surface Patches"          %
%                                                                         %
%*************************************************************************%
% This library is free software; you can redistribute it and/or modify	  %
% it under the terms of the GNU Lesser General Public License as published%
% by the Free Software Foundation; either version 2.1 of the License, or  %
% (at your option) any later version.                                     %
%                                                                         %
% This library is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU Lesser General Public License for more details.        	  %
% You should have received a copy of the GNU Lesser General Public License%
% along with this library; if not, write to the Free Software Foundation, %
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA       %
%*************************************************************************%
% Author of this file: Mahmoud S. Shaqfa
% The Pareto-like Sequential Sampling Algorithm (PSS)
% This is a simple implementation of the algorithm
% The main code was implemented in C++14 and Python.
% If you use this code please cite this paper:
%       Shaqfa, M., Beyer, K. Pareto-like sequential sampling heuristic for 
%       global optimisation. Soft Comput (2021).
%       https://doi.org/10.1007/s00500-021-05853-8

function results_struct = PSS_optimizer(fun, LB, UB, pop_size, iterations, acceptance_rate, plot_results, silent_mode)
tic; format long;
if nargin == 7
    silent_mode =  false;
end
%% The objective function parameters
dim = length(LB);
exact = 0.;
%% Initialize and pre-allocation of memory containers
% LB = ones(1, dim) * LB_val; % Lower bounds vector
% UB = ones(1, dim) * UB_val; % Upper bounds vector
funval_his = zeros(iterations + 1, 1); % The history of all iterations
funval = zeros(pop_size, 1); % Function evaluations per iteration for all the population memebers
limits_his = zeros(iterations + 1, dim * 2); % The histories of LB and UB

%% Initialize the population - Fast vectorized initial sampling
LB_mat = repmat(LB, pop_size, 1);
UB_mat = repmat(UB, pop_size, 1);
random_pop = rand(pop_size, dim);
population = LB_mat + random_pop .* (UB_mat - LB_mat);

clear LB_mat UB_mat % we don't need them after this stage (clear memory)

% Evaluate the inital population
for i = 1: pop_size
    funval(i) = fun(population(i, :));
end

% Find the current best solution
[~, best_index] = min(funval);
funval_his(1) = funval(best_index);
best_sol = population(best_index, :);
new_solution = true;
%% The core of the algorithm (the critical loop)
for i = 2: iterations + 1
    random_pop = rand(pop_size, dim);
    for j = 1:pop_size
       counter = 1;
       for k = 1: dim
           % Change the prominent domain (a self-adaptive mechanism)
           if new_solution
              deviation = 0.5 * (1 - acceptance_rate) * (1 - (i/iterations)) ...
                  * (UB(k) - LB(k));
           end
           reduced_LB = max([best_sol(k) - deviation, LB(k)]);
           reduced_UB = min([reduced_LB + 2 * deviation, UB(k)]);
           
           % Hold the history of the prominent domain
           limits_his(i, counter) = reduced_LB;
           counter = counter + 1;
           limits_his(i, counter) = reduced_UB;
           counter = counter + 1;
           
           if rand <= acceptance_rate
               % Sample from the prominent domain
               population(j, k) = reduced_LB + random_pop(j, k) ...
                   * (reduced_UB - reduced_LB);
           else
               % Sample from the overall domain
               population(j, k) = LB(k) + random_pop(j, k) ...
                   * (UB(k) - LB(k));
           end
       end
       % Evaluate the solution vector
       funval(j) = fun(population(j, :));    
    end
    % Find the current best value
    if min(funval) <= funval_his(i - 1)
        [~, best_index] = min(funval);
        best_sol = population(best_index, :);
        new_solution = true;
        funval_his(i) = min(funval);
    else
        new_solution = false;
        funval_his(i) = funval_his(i - 1);
    end
    if ~silent_mode
        fprintf('\nIteration number %2.0f and best fitness is: %2.10f', i, funval_his(i))
    end
end
time = toc;
if ~silent_mode
    fprintf('\n\nThe elapsed time %4.5f min', time / 60)
    fprintf('\n\nThe best solution is: \n')
    disp(best_sol)
end
%% Plot the results
if ~silent_mode
    disp('Exporting the figures, please wait ...')
end
if plot_results
    % plot convergence
    f1 = figure;
    plot(1:length(funval_his), funval_his, 'LineWidth', 2)
    xlim([2, length(funval_his)])
    xlabel('Iterations')
    ylabel('Area distortion')
    title('Convergence curve')
    print(f1, 'convergence', '-r400', '-dpdf')
    set(gca,'fontname','Amiri')  % Set it to times
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    % plot the prominent domain history
    f2 = figure('visible','off');
    set(f2, 'Position', [0 0 200 20*dim])
    counter = 1;
    for kk = 1:dim
        subplot(ceil(dim/2), 2, kk)
        title(num2str(kk))
        plot(1:length(limits_his(:,1)), limits_his(:, counter), '-r')
        counter = counter + 1;
        hold on
        plot(1:length(limits_his(:,1)), limits_his(:, counter), '-b')
        counter = counter + 1;
        hold on
        plot(1:length(limits_his(:,1)), ...
            ones(1, length(limits_his(:,1)))*exact, '--m')
        xlim([2, length(funval_his)])
        hold off
    end
    saveas(f2,'prominent_domain','pdf')
    print(f2, 'prominent_domain', '-r400', '-dpdf','-fillpage')
end
if ~silent_mode
    disp('Done')
end
%% Saving the final results and vectors
results_struct.time = time;
results_struct.evaluations_history = funval_his;
results_struct.LB_UB_history = limits_his;
results_struct.best_solution = best_sol;
end