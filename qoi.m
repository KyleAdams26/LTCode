function QOI = qoi(p, debug)
%Author: Dr. Jaimit Parikh
    arguments
        p; % parameters
        debug = false; % flag to plot for debugging purpose
    end
    
    t0 = 0; tfinal = 30; % simulation time in days
    IC = getInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];
    % the state variables

    %simulates the model 
    [T,Y] = ode15s(@(t, y)odefun(t, y, p),...
        [t0 tfinal], IC);
    %we chose the L(T) as our QOI   
    QOI = Y(end,1);
    % results
    if debug
    plotQOI(T, Y);
    end
end
%this function is a plot used for debugging
function plotQOI(T, Y, fname)
    f = figure('DefaultAxesFontSize',12,...
        'Position', [20 20 1800 900]);
    sp = ["L", "A", "Tc", "Th", "Tr", "I"];
    for idx = 1:size(Y, 2)
        subplot(3, 5, idx);
        plot(T, Y(:, idx), 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'matlab');
        xlabel('Time, days');
        ylabel(sp(idx));
        xlim([0, 30]);
    end
    legend();
    exportgraphics(f, fname, 'resolution', 300);
end
