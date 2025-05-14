function modelSimulationsWithParamAdj
    
    %initializing colors
    L_color = [166/255, 107/255, 97/255];

    %loading parameters
    param = parameters();

    %solve model regularly
    t0 = 0; tfinal = 30; % simulation time in days
    IC = getInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];
    tspan = [t0 tfinal];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12) ;

    % the state variables

    %simulates the model 
    [T,Y] = ode45(@(t, y)odefun(t, y, param),...
        [t0 tfinal], IC, options);
    L = Y(:, 1);
    %choose parameters to show their influence on QOI
    [T1, Y1, T2, Y2] = solveWithModifiedParam('dL', param, tspan); L1 = Y1(:,1); L2 = Y2(:,1);
    [T3, Y3, T4, Y4] = solveWithModifiedParam('aCL', param, tspan); L3 = Y3(:,1); L4 = Y4(:,1);
    [T5, Y5, T6, Y6] = solveWithModifiedParam('gC', param, tspan); L5 = Y5(:,1); L6 = Y6(:,1);
    
    %plot simulations
    figure;
    tiledlayout(2, 3)

    addTile('dL', T1, L1, T2, L2, L_color, 30);
    addTile('aCL', T3, L3, T4, L4, L_color, 30);
    addTile('gC', T5, L5, T6, L6, L_color, 30);





   function addTile(paramName, x_vals1_5, nominal_param1_5, x_vals_half, nominal_param_half, graph_color, endTime)
    nexttile;
    hold on;
    plot(T, L, 'Color', graph_color, 'LineWidth', 3);
    plot(x_vals1_5, nominal_param1_5, 'Color', graph_color, 'LineWidth', 3, 'LineStyle', "--");
    plot(x_vals_half, nominal_param_half, 'Color', graph_color, 'LineWidth', 3, 'LineStyle', ":");
    
    %Axis settings
    xlim([0 endTime]);
    xticks(0:10:30);
    ylim([0.4e11 2e11]);
    xlabel('Time (Days)');
    ylabel('L (cells/\mu L)');
    title('Liver hepatocytes (L)');
    
    %Font size adjustments
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    
    %Legend
    legend({sprintf('%s', paramName),...
        sprintf('1.5*%s', paramName),...
        sprintf('0.5*%s', paramName)}, 'Location', 'southwest', 'FontSize', 10, 'Box', 'off');
    hold off;
end

function [T1, Y1, T2, Y2] = solveWithModifiedParam(paramName, p, tspan)
    % Solve ODE with modified parameter values (1.5x and 0.5x of given param).
    % Inputs:
    %   paramName - String of the parameter to modify
    %   p - Struct containing model parameters
    %   tspan - Time span for ODE solver
    % Outputs:
    %   T1, Y1 - Solution for 1.5 * given parameter
    %   T2, Y2 - Solution for 0.5 * given parameter

    %get initial conditions
    IC = getInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];
    %options
    options = odeset('RelTol',1e-12,'AbsTol',1e-12) ;

    %modify parameter and solve with 1.5x value
    p1 = p;
    p1.(paramName) = 1.5 * p.(paramName);
    [T1, Y1] = ode45(@(t, y) odefun(t, y, p1), tspan, IC, options);

    %modify parameter and solve with 0.5x value
    p2 = p;
    p2.(paramName) = 0.5 * p.(paramName);
    [T2, Y2] = ode45(@(t, y) odefun(t, y, p2), tspan, IC);
end


end 
