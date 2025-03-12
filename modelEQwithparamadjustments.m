function modelEQwithparamadjustments
    
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
    gC = 2.08; %12 INFLUENTIAL
    KC = 598; %13
    aIC = 2; %14 INFLUENTIAL
    bIC = .178; %15
    lH = 0.0000063; %16
    gH = 1.51; %17
    KH = 422; %18
    aIH = 2; %19
    bIH = 0.178; %20
    lR = .000000667; %21
    dI = 166; %22 INFLUENTIAL
    aHC = 1; %23
    bHC = 35; %24 
    dC = 0.585*1; %25 INFLUENTIAL
    
    aAH = 0.0000261; %26
    bAH = 4; %27
    aRA = 0.4; %28
    bRA = 20; %29
    aIRA = 2; %30
    bIRA = 0.356; %31
    dH = 0.333; %32

    dL = 0.005; %33 INFLUENTIAL
    aCL = 10; %34 INFLUENTIAL
    bCL = 200; %35 INFLUENTIAL

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

    %solving with different dCs
    [Tf, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL, dI, aIC, gC), tspan, x0, options);
    ogLF = Xf(:,6);
    [dCTf2, Xf] = ode45(@(t, y) Forward(t, y, dC*2, dL, aCL, bCL, dI, aIC, gC), tspan, x0, options);
    dC_LF2 = Xf(:,6);
    [dCTf_half, Xf] = ode45(@(t, y) Forward(t, y, dC*0.5, dL, aCL, bCL, dI, aIC, gC), tspan, x0, options);
    dC_LF_half = Xf(:,6);

    %solving with different dLs
    [dLTf2, Xf] = ode45(@(t, y) Forward(t, y, dC, dL*2, aCL, bCL, dI, aIC, gC), tspan, x0, options);
    dL_LF2 = Xf(:,6);
    [dLTf_half, Xf] = ode45(@(t, y) Forward(t, y, dC, dL*0.5, aCL, bCL, dI, aIC, gC), tspan, x0, options);
    dL_LF_half = Xf(:,6);

     %solving with different aCLs
    [aCLTf2, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL*2, bCL, dI, aIC, gC), tspan, x0, options);
    aCL_LF2 = Xf(:,6);
    [aCLTf_half, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL*0.5, bCL, dI, aIC, gC), tspan, x0, options);
    aCL_LF_half = Xf(:,6);

    %solving with different bCLs
    [bCLTf2, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL*2, dI, aIC, gC), tspan, x0, options);
    bCL_LF2 = Xf(:,6);
    [bCLTf_half, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL*0.5, dI, aIC, gC), tspan, x0, options);
    bCL_LF_half = Xf(:,6);

    %solving with different dIs
    [dITf2, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL, dI*2, aIC, gC), tspan, x0, options);
    dI_LF2 = Xf(:,6);
    [dITf_half, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL, dI*0.5, aIC, gC), tspan, x0, options);
    dI_LF_half = Xf(:,6);

    %solving with different aICs
    [aICTf2, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL, dI, aIC*2, gC), tspan, x0, options);
    aIC_LF2 = Xf(:,6);
    [aICTf_half, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL, dI, aIC*0.5, gC), tspan, x0, options);
    aIC_LF_half = Xf(:,6);

    %solving with different gCs
    [gCTf2, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL, dI, aIC, gC*2), tspan, x0, options);
    gC_LF2 = Xf(:,6);
    [gCTf_half, Xf] = ode45(@(t, y) Forward(t, y, dC, dL, aCL, bCL, dI, aIC, gC*0.5), tspan, x0, options);
    gC_LF_half = Xf(:,6);

    %all for /  all against
    [all_for_Tf, Xf] = ode45(@(t, y) Forward(t, y, dC*2, dL*0.5, aCL*0.5, bCL*2, dI*2, aIC*0.5, gC*0.5), tspan, x0, options);
    all_for_LF = Xf(:,6);
    [all_against_Tf, Xf] = ode45(@(t, y) Forward(t, y, dC*0.5, dL*2, aCL*2, bCL*0.5, dI*0.5, aIC*2, gC*2), tspan, x0, options);
    all_against_LF = Xf(:,6);

    %% Plots Together %%

    tiledlayout(3,3)

    nexttile
    %Plot for dC
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(dCTf2, dC_LF2, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(dCTf_half, dC_LF_half, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'dC', '2*dC', '0.5*dC'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')
    hold off;

    %dL plot
    nexttile
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(dLTf2, dL_LF2, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(dLTf_half, dL_LF_half, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'dL', '2*dL', '0.5*dL'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')

    hold off;

    %aCL plot
    nexttile
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(aCLTf2, aCL_LF2, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(aCLTf_half, aCL_LF_half, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'aCL', '2*aCL', '0.5*aCL'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')
    hold off;
    
    %bCL
    nexttile
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(bCLTf2, bCL_LF2, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(bCLTf_half, bCL_LF_half, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'bCL', '2*bCL', '0.5*bCL'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')
    hold off;

    %dI
    nexttile
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(dITf2, dI_LF2, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(dITf_half, dI_LF_half, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'dI', '2*dI', '0.5*dI'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')
    hold off;

    %aIC
    nexttile
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(aICTf2, aIC_LF2, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(aICTf_half, aIC_LF_half, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'aIC', '2*aIC', '0.5*aIC'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')
    hold off;
    
    %gC
    nexttile
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(gCTf2, gC_LF2, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(gCTf_half, gC_LF_half, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'gC', '2*gC', '0.5*gC'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')
    hold off;
    
    %all at once
    nexttile
    hold on;
    plot(Tf, ogLF, 'Color', L_color, 'LineWidth',3)
    plot(all_for_Tf, all_for_LF, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "--")
    plot(all_against_Tf, all_against_LF, 'Color', L_color, 'Linewidth', 3, 'Linestyle', "-.")
    xlim([0 endTime])
    xticks(0: 10 : 30)
    xlabel('Time (Days)')
    ylabel('L (cells/\mu L)')
    title('Liver hepatocytes (L) over time')
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    legend({'Nominal', 'All for', 'All against'}, 'Location', 'southwest', 'FontSize', 7, 'Box', 'off')
    hold off;
    


    function dy = Forward(~,y, dCparam, dLparam, aCLparam, bCLparam, dIparam, aICparam, gCparam)
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
    pPath = gCparam*Tc*(1 - Tc/KC);
    qPath = aICparam*I/(bIC + I);
    sPath = gH*Th*(1 - Th/KH);
    rPath = aIH*I/(bIH + I);
    wPath = dIparam*I;
    lPath = aHC*Th/(bHC + Th);
    tPath = dCparam*Tc;
    gPath = aAH*A/(bAH + A);
    jPath = aRA*Tr/(bRA + Tr);
    xPath = 1 + aIRA*I/(bIRA + I);
    yPath = dH*Th;
    ePath = aCLparam*Tc/(bCLparam + Tc);


    % Dynamics 
    dy(1)  = lL*(dLparam)*L*(1+ePath) - iPath;  %dAdt
    dy(2)  = vPath - zPath*uPath; %dTrdt
    dy(3)  = nPath + oPath -lC*pPath*qPath -lH*sPath*rPath - lR*zPath*uPath - wPath;
    dy(4)  = lPath + pPath*qPath - tPath; %dTcdt
    dy(5)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;
    dy(6) = -dLparam*L*ePath; %dLdt
        

    end

end 
