function QOI = calculateQOI(p)
%Original code by Jaimit Parikh; modified by Kyle Adams
%This file saves the QOI values after running the ODE model
    arguments
        p; % parameters
    end
    
    % ###Step 1
    t0 = 0; tfinal = 30; % simulation time in days
    IC = setInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];

    %simulate the model 
    [~,Y] = ode15s(@(t, y)odefun(t, y, p),...
        [t0 tfinal], IC); 
    %we chose the L(T) as our QOI   
    QOI = Y(end,1); % ###Step 2
    % results
end
