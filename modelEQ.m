function modelEQ
    
    endTime = 30; %initializing QOI time
    %initializing colors
    A_color = [250/255, 129/255, 113/255]; 
    Tr_color = [139/255, 235/255, 229/255];
    I_color = [251/255, 186/255, 27/255];
    Tc_color = [252/255, 132/255, 217/255];
    Th_color = [176/255, 101/255, 243/255];
    L_color = [166/255, 107/255, 97/255];

    % =================================================================== %
    % ------ Initialize constant parameters ------ %
    global lL dA ...
            sR dR aIR bIR ...
            aCI bCI aHI bHI lC gC KC aIC bIC lH gH KH bIH lR dI ...
            aHC bHC dC ...
            aAH bAH aRA bRA aIRA bIRA dH ...
            dL aCL bCL ... %35 parameters

    % ------ Initialize parameters  ------ %
    lL = 0.00000000452; %1
    dA = 0.0833; %2

    sR = 0.107; %3
    dR = 0.0658; %4
    aIR = 0.625; %5
    bIR = 0.00833; %6
    
    aCI = 0.36; %7
    bCI = 352; %8
    aHI = 70.7; %9
    bHI = 99.7; %10
    lC = 0.001; %11
    gC = 2.08; %12
    KC = 598; %13
    aIC = 2; %14
    bIC = .178; %15
    lH = 0.0000063; %16
    gH = 1.51; %17
    KH = 422; %18
    aIH = 2; %19
    bIH = 0.178; %20
    lR = .000000667; %21
    dI = 166; %22
    aHC = 1; %23
    bHC = 35; %24 
    dC = 0.585; %25
    
    aAH = 0.0000261; %26
    bAH = 4; %27
    aRA = 0.4; %28
    bRA = 20; %29
    aIRA = 2; %30
    bIRA = 0.356; %31
    dH = 0.333; %32

    dL = 0.005; %33
    aCL = 10; %34
    bCL = 200; %35

    % Initial values for variables 
    A0 = 8; % Alloantigen
    Tr0 = 1.31; % Treg cells
    I0 = 0.0113; % IL-2
    Tc0 = 100; % Tc cells  
    Th0 = 70; % Helper T cells
    L0 = 2*(10^11); % Liver cells


    % The integration settings defined by options
    options = odeset('RelTol',1e-12,'AbsTol',1e-12) ;
    x0 = [A0, Tr0, I0, Tc0, Th0, L0] ;
    tspan = [0 endTime] ; 
    [Tf, Xf] = ode45(@Forward, tspan, x0, options); %solving system
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
    xlim([0 endTime])
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
    xlim([0 endTime]) 
    xlabel('Time (Days)')
    ylabel('A (cells/\mu L)')
    title('APCs (A) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Top right plot
    nexttile
    plot(Tf,TcF, 'Color', Tc_color, 'LineWidth', 3)
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('T_{C} (cells/\mu L)')
    title('Cytotoxic T cells (T_{C}) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Bottom left plot
    nexttile
    plot(Tf,ThF, 'Color', Th_color,'LineWidth',3 )
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('T_{H} (cells/\mu L)')
    title('Helper T cells (T_{H}) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Bottom middle plot
    nexttile
    plot(Tf,TrF,'Color', Tr_color, 'LineWidth',3)
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('T_{R} (cells/\mu L)')
    title('Regulatory T cells (T_{R}) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;

    % Bottom right plot
    nexttile
    plot(Tf,IF,'Color', I_color,'LineWidth',3)
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('I (ng/\mu L)')
    title('IL-2 (I) over time')
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
   

    %the following function represents the system of ODEs
    function dy = Forward(~,y)
        dy = zeros(6,1);
        % States 
        A = y(1) ;
        Tr = y(2) ;
        I = y(3) ;
        Tc = y(4) ;
        Th = y(5) ;
        L = y(6) ;
   
        %% -- The rate change of the 6 populations (Dynamics of the system) -- %%
        % Paths for the dynamics %
   
    iPath = dA*A;
    vPath = sR;
    zPath = dR*Tr;
    uPath = 1 - aIR*I/(bIR + I);
    nPath = aCI*Tc/(bCI +Tc);
    oPath = aHI*Th/(bHI + Th);
    pPath = gC*Tc*(1 - Tc/KC);
    qPath = aIC*I/(bIC + I);
    sPath = gH*Th*(1 - Th/KH);
    rPath = aIH*I/(bIH + I);
    wPath = dI*I;
    lPath = aHC*Th/(bHC + Th);
    tPath = dC*Tc;
    gPath = aAH*A/(bAH + A);
    jPath = aRA*Tr/(bRA + Tr);
    xPath = 1 + aIRA*I/(bIRA + I);
    yPath = dH*Th;
    ePath = aCL*Tc/(bCL + Tc);


    % Dynamics 
    dy(1)  = lL*(dL)*L*(1+ePath) - iPath;  %dAdt
    dy(2)  = vPath - zPath*uPath; %dTrdt
    dy(3)  = nPath + oPath -lC*pPath*qPath -lH*sPath*rPath - lR*zPath*uPath - wPath;
    dy(4)  = lPath + pPath*qPath - tPath; %dTcdt
    dy(5)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;
    dy(6) = -dL*L*ePath; %dLdt
        

    end
end 
