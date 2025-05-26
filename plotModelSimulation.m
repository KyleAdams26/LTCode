%%adjusted by Kyle Adams, originally Mahya Aghaee

function plotModelSimulation
    % ###Step 1
    % define the colors of the state variables
    L_color  = [166/255, 107/255,  97/255]; % brown
    A_color  = [250/255, 129/255, 113/255]; % salmon
    Th_color = [176/255, 101/255, 243/255]; % purple
    Tc_color = [252/255, 132/255, 217/255]; % pink
    Tr_color = [139/255, 235/255, 229/255]; % turquoise
    I_color  = [251/255, 186/255,  27/255]; % goldenrod
    line_width = 4;

    %load the parameter names and values
    p = setParameters();

    % ###Step 2
    % set simulation timespan and load the initial conditions
    t0 = 0; tfinal = 30; % initial and final simulation times in days
    IC = setInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];

    % ###Step 3
    % choose integration settings using "options"
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    % simulate the model 
    [Tf,Xf] = ode45(@(t, y)odefun(t, y, p),...
        [t0 tfinal], IC, options);
    
    % ###Step 4
    % store the solutions for each variable for plotting
    LF   = Xf(:,1);
    AF   = Xf(:,2);
    ThF  = Xf(:,3);
    TcF  = Xf(:,4);
    TrF  = Xf(:,5);
    IF   = Xf(:,6);
   
    % ###Step 5
    %% Create plots of the variables in a 2x3 grid %%
    % each plot has the same aesthetic and similar labelings
    
    tiledlayout(2,3)

    % Top left plot
    nexttile
    plot(Tf, LF, 'Color', L_color, 'LineWidth',line_width)
    xlim([0 tfinal])
    xlabel('Time (days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    formatAxes(ax);

    % Top middle plot
    nexttile
    plot(Tf,AF,'Color', A_color,'LineWidth', line_width)
    title('Cell Populations')
    xlim([0 tfinal]) 
    xlabel('Time (days)')
    ylabel('A (cells/\mu L)')
    title('APCs (A) over time')
    ax = gca;
    formatAxes(ax);

    % Top right plot
    nexttile
    plot(Tf,TcF, 'Color', Tc_color, 'LineWidth', line_width)
    xlim([0 tfinal])
    xlabel('Time (days)')
    ylabel('T_{C} (cells/\mu L)')
    title('Cytotoxic T cells (T_{C}) over time')
    ax = gca;
    formatAxes(ax);

    % Bottom left plot
    nexttile
    plot(Tf,ThF, 'Color', Th_color,'LineWidth',line_width)
    xlim([0 tfinal])
    xlabel('Time (days)')
    ylabel('T_{H} (cells/\mu L)')
    title('Helper T cells (T_{H}) over time')
    ax = gca;
    formatAxes(ax);

    % Bottom middle plot
    nexttile
    plot(Tf,TrF,'Color', Tr_color, 'LineWidth',line_width)
    xlim([0 tfinal])
    xlabel('Time (days)')
    ylabel('T_{R} (cells/\mu L)')
    title('Regulatory T cells (T_{R}) over time')
    ax = gca;
    formatAxes(ax);

    % Bottom right plot
    nexttile
    plot(Tf,IF,'Color', I_color,'LineWidth',line_width)
    xlim([0 tfinal])
    xlabel('Time (days)')
    ylabel('I (ng/\mu L)')
    title('IL-2 (I) over time')
    ax = gca;
    formatAxes(ax);
  
    %helper function for formatting axes
    function formatAxes(ax)
        ax.Title.FontSize = 15;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
    end

end 
