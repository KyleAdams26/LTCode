function modelEQ
    
    endTime = 30; %initializing QOI time
    %initializing colors
    A_color = [254/255, 78/255, 76/255]; 
    Tr_color = [139/255, 235/255, 229/255];
    I_color = [251/255, 186/255, 27/255];
    Tc_color = [252/255, 132/255, 217/255];
    Th_color = [176/255, 101/255, 243/255];
    L_color = [166/255, 107/255, 97/255];

    % =================================================================== %
    % ------ Initialize constant parameters ------ %
    global aLA dA ...
            sR dR aIR bIR ...
            aCI bCI aHI bHI lC gC KC bIC lH gH KH bIH lR dI ...
            aHC bHC dC ...
            aAH bAH aRA bRA aIRA bIRA dH ...
            dL aCL bCL ... %33 parameters

    % ------ Initialize parameters  ------ %
    aLA = 1; %1
    dA = 0.083333; %2

    sR = 0.107; %3
    dR = 0.0658; %4
    aIR = 0.625; %5
    bIR = 0.00833333; %6
    
    aCI = 0.36; %7
    bCI = 350; %8
    aHI = 70.7; %9
    bHI = 100; %10
    lC = 0.005; %11
    gC = 2.079; %12
    KC = .2*2991; %13
    bIC = .178; %14
    lH = 0.000005999988; %15
    gH = 1.512; %16
    KH = .2*2112; %17
    bIH = 0.178; %18
    lR = .00000333; %19
    dI = 166.355; %20
    
    aHC = 1; %21
    bHC = 35; %22 
    dC = 0.5853658537; %23
    
    aAH = 0.0000261; %24
    bAH = 4; %25
    aRA = 0.4; %26
    bRA = 20; %27
    aIRA = 2; %28
    bIRA = 0.356; %29
    dH = 0.3333; %30

    dL = 1/200; %31
    aCL = 10; %32
    bCL = 200; %33

    % Initial values for variables 
    A0 = 8; % Alloantigen
    Tr0 = 1.307; % Treg cells
    I0 = 0.00565; % IL-2
    Tc0 = 100; % Tc cells  
    Th0 = 70; % Helper T cells
    L0 = 2*(10^9); % Liver cells


    % The integration settings defined by options
    options = odeset('RelTol',1e-12,'AbsTol',1e-12) ;
    x0 = [A0, Tr0, I0, Tc0, Th0, L0] ;
    tspan = [0 endTime] ; 
    [Tf, Xf] = ode45(@Forward, tspan, x0, options);
    AF  = Xf(:,1);
    TrF = Xf(:,2);
    IF  = Xf(:,3);
    TcF = Xf(:,4);
    ThF = Xf(:,5);
    LF = Xf(:,6);
   
    %% Plots Together %%

    tiledlayout(2,3)

    % Top left plot
    nexttile
    plot(Tf, LF, 'Color', L_color, 'LineWidth',3)
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('Liver hepatocytes')
    title('Liver hepatocytes (L) over time')

    % Top middle plot
    nexttile
    plot(Tf,AF,'Color', A_color,'LineWidth', 3)
    title('Cell Populations')
    xlim([0 endTime]) 
    xlabel('Time (Days)')
    ylabel('APCs')
    title('APCs (A) over time')

    % Top right plot
    nexttile
    plot(Tf,TcF, 'Color', Tc_color, 'LineWidth', 3)
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('Cytotoxic T cells')
    title('Cytotoxic T cells (Tc) over time')

    % Bottom left plot
    nexttile
    plot(Tf,ThF, 'Color', Th_color,'LineWidth',3 )
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('Helper T cells')
    title('Helper T cells (Th) over time')

    % Bottom middle plot
    nexttile
    plot(Tf,TrF,'Color', Tr_color, 'LineWidth',3)
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('Regulatory T cells')
    title('Regulatory T cells (Tregs) over time')

    % Bottom right plot
    nexttile
    plot(Tf,IF,'Color', I_color,'LineWidth',3)
    xlim([0 endTime])
    xlabel('Time (Days)')
    ylabel('IL-2')
    title('IL-2 (I) over time')
   

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
    qPath = I/(bIC + I);
    sPath = gH*Th*(1 - Th/KH);
    rPath = I/(bIH + I);
    wPath = dI*I;
    lPath = aHC*Th/(bHC + Th);
    tPath = dC*Tc;
    gPath = aAH*A/(bAH + A);
    jPath = aRA*Tr/(bRA + Tr);
    xPath = 1 + aIRA*I/(bIRA + I);
    yPath = dH*Th;
    ePath = aCL*Tc/(bCL + Tc);


    % Dynamics 
    dy(1)  = aLA*(dL)*L*(1+ePath) - iPath ; %dAdt
    dy(2)  = vPath - zPath*uPath; %dTrdt
    dy(3)  = nPath + oPath -lC*pPath*qPath -lH*sPath*rPath - lR*zPath*uPath - wPath;
    dy(4)  = lPath + pPath*qPath - tPath; %dTcdt
    dy(5)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;
    dy(6) = -dL*L*ePath; %dLdt
        

    end
end 