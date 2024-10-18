function dcdt = odefun(~, c, p)
    % odefun - system of equations of the Liver Transplant model
    % Inputs:
    %   ~ - time
    %   c - state variables
    %   p - parameter structure

    % Outputs:
    % dcdt - Derivatives of the state variables
    % The credit of this code is for Dr. Jaimit Parikh
    % We modified his code to be suitable for our model. 
    % This code was originally modified by Mahya Aghaee.

    % # c : levels of 6 biological factors at t. 
    A = c(1); % Alloantigen
    Tr = c(2); % Treg cells
    I = c(3); % IL-2
    Tc = c(4); % Tc cells  
    Th = c(5); % Helper T cells
    L = c(6); % Liver cells

    %     # t : time (int)
    %     # x : 33 parameter values

    %     #   prepare parameter values
    aLA = p.aLA; %1
    bLA = p.bLA; %2
    dA = p.dA; %3

    sR = p.sR; %4
    dR = p.dR; %5
    aIR = p.aIR; %6
    bIR = p.bIR; %7
    
    aCI = p.aCI; %8
    bCI = p.bCI; %9
    aHI = p.aHI; %10
    bHI = p.bHI; %11
    lC = p.lC; %12
    gC = p.gC; %13
    KC = p.KC; %14
    bIC = p.bIC; %15
    lH = p.lH; %16
    gH = p.gH; %17
    KH = p.KH; %18
    bIH = p.bIH; %19
    lR = p.lR; %20
    dI = p.dI; %21
    
    aHC = p.aHC; %22
    bHC = p.bHC; %23
    dC = p.dC; %24
    
    aAH = p.aAH; %25
    bAH = p.bAH; %26
    aRA = p.aRA; %27
    bRA = p.bRA; %28
    aIRA = p.aIRA; %29
    bIRA = p.bIRA; %30
    dH = p.dH; %31

    dL = p.dL; %32
    aCL = p.aCL; %33
    bCL = p.bCL; %34




        %% -- The rate change of the 6 populations (Dynamics of the system) -- %%
        % Paths for the dynamics %
   
    bPath = aLA*L/(bLA + L);
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
    dy(1)  = bPath - iPath ;
    dy(2)  = vPath - zPath*uPath;
    dy(3)  = nPath + oPath -lC*pPath*qPath -lH*sPath*rPath - lR*zPath*uPath - wPath;
    dy(4)  = lPath + pPath*qPath - tPath;
    dy(5)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;
    dy(6) = -dL*L*ePath;
        


    dcdt = [dy(1), dy(2), dy(3), dy(4), dy(5), dy(6)]';
end