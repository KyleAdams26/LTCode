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
    lL = p.lL; %1
    dA = p.dA; %2

    sR = p.sR; %3
    dR = p.dR; %4
    aIR = p.aIR; %5
    bIR = p.bIR; %6
    
    aCI = p.aCI; %7
    bCI = p.bCI; %8
    aHI = p.aHI; %9
    bHI = p.bHI; %10
    lC = p.lC; %11
    gC = p.gC; %12
    KC = p.KC; %13
    aIC = p.aIC; %14
    bIC = p.bIC; %15
    lH = p.lH; %16
    gH = p.gH; %17
    KH = p.KH; %18
    aIH = p.aIH; %19
    bIH = p.bIH; %20
    lR = p.lR; %21
    dI = p.dI; %22
    
    aHC = p.aHC; %23
    bHC = p.bHC; %24
    dC = p.dC; %25
    
    aAH = p.aAH; %26
    bAH = p.bAH; %27
    aRA = p.aRA; %28
    bRA = p.bRA; %29
    aIRA = p.aIRA; %30
    bIRA = p.bIRA; %31
    dH = p.dH; %32

    dL = p.dL; %33
    aCL = p.aCL; %34
    bCL = p.bCL; %35




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
    dy(1)  = lL*(dL)*L*(1+ePath) - iPath ; %dAdt
    dy(2)  = vPath - zPath*uPath; %dTrdt
    dy(3)  = nPath + oPath -lC*pPath*qPath -lH*sPath*rPath - lR*zPath*uPath - wPath;
    dy(4)  = lPath + pPath*qPath - tPath; %dTcdt
    dy(5)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;
    dy(6) = -dL*L*ePath; %dLdt
        


    dcdt = [dy(1), dy(2), dy(3), dy(4), dy(5), dy(6)]';
end
