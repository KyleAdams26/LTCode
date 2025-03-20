function modelEQ
    %initializing colors
    A_color = [250/255, 129/255, 113/255]; 
    Tr_color = [139/255, 235/255, 229/255];
    I_color = [251/255, 186/255, 27/255];
    Tc_color = [252/255, 132/255, 217/255];
    Th_color = [176/255, 101/255, 243/255];
    L_color = [166/255, 107/255, 97/255];

    %load in parameters
    p = parameters();

    %choosing timespan and initial conditions
    t0 = 0; tfinal = 30; % simulation time in days
    IC = getInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];

   
    % The integration settings defined by options
    options = odeset('RelTol',1e-12,'AbsTol',1e-12) ;
    %simulates the model 
    [Tf,Xf] = ode45(@(t, y)odefun(t, y, p),...
        [t0 tfinal], IC, options);
    %the next lines store the solutions for each variable of the system for plotting
    AF  = Xf(:,1);
    TrF = Xf(:,2);
    IF  = Xf(:,3);
    TcF = Xf(:,4);
    ThF = Xf(:,5);
    LF = Xf(:,6);
   
    %% Plots Together %%
    % each plot has the same aesthetic and similar labelings
    
    tiledlayout(2,3)

    % Top left plot
    nexttile
    plot(Tf, LF, 'Color', L_color, 'LineWidth',3)
    xlim([0 tfinal])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Top middle plot
    nexttile
    plot(Tf,AF,'Color', A_color,'LineWidth', 3)
    title('Cell Populations')
    xlim([0 tfinal]) 
    xlabel('Time (Days)')
    ylabel('A (cells/\mu L)')
    title('APCs (A) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Top right plot
    nexttile
    plot(Tf,TcF, 'Color', Tc_color, 'LineWidth', 3)
    xlim([0 tfinal])
    xlabel('Time (Days)')
    ylabel('T_{C} (cells/\mu L)')
    title('Cytotoxic T cells (T_{C}) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Bottom left plot
    nexttile
    plot(Tf,ThF, 'Color', Th_color,'LineWidth',3 )
    xlim([0 tfinal])
    xlabel('Time (Days)')
    ylabel('T_{H} (cells/\mu L)')
    title('Helper T cells (T_{H}) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Bottom middle plot
    nexttile
    plot(Tf,TrF,'Color', Tr_color, 'LineWidth',3)
    xlim([0 tfinal])
    xlabel('Time (Days)')
    ylabel('T_{R} (cells/\mu L)')
    title('Regulatory T cells (T_{R}) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Bottom right plot
    nexttile
    plot(Tf,IF,'Color', I_color,'LineWidth',3)
    xlim([0 tfinal])
    xlabel('Time (Days)')
    ylabel('I (ng/\mu L)')
    title('IL-2 (I) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
  
end 