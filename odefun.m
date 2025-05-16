function dcdt = odefun(~, c, p) %QC'd
    % odefun - system of equations of the Liver Transplant model
    % Inputs:
    %   ~ - time
    %   c - state variables
    %   p - parameter structure
    % Outputs:
    % dcdt - Derivatives of the state variables
    % This code is based on code originally written by Jaimit Parikh
    % in 2023 and modified by Mahya Aghaee in 2024 
    % Kyle Adams modified this code in 2025 to be suitable for the Liver Transplant
    % Project.
 
    % # c : levels of 6 biological factors at t. 
    L  = c(1); % Liver cells
    A  = c(2); % APCs
    Th = c(3); % Helper T cells
    Tc = c(4); % Tc cells  
    Tr = c(5); % T reg cells
    I  = c(6); % IL-2

    %     # t : time (int)
    %     # x : 35 parameter values

        %% -- The rate of change of the 6 populations (Dynamics of the system) -- %%
        % Pathways for the dynamics %
   
    iPath = p.dA*A;
    vPath = p.sR;
    zPath = p.dR*Tr;
    uPath = 1 - p.aIR*I/(p.bIR + I);
    nPath = p.aCI*Tc/(p.bCI +Tc);
    oPath = p.aHI*Th/(p.bHI + Th);
    pPath = p.gC*Tc*(1 - Tc/p.KC);
    qPath = p.aIC*I/(p.bIC + I);
    sPath = p.gH*Th*(1 - Th/p.KH);
    rPath = p.aIH*I/(p.bIH + I);
    wPath = p.dI*I;
    lPath = p.aHC*Th/(p.bHC + Th);
    tPath = p.dC*Tc;
    gPath = p.aAH*A/(p.bAH + A);
    jPath = p.aRA*Tr/(p.bRA + Tr);
    xPath = 1 + p.aIRA*I/(p.bIRA + I);
    yPath = p.dH*Th;
    ePath = p.aCL*Tc/(p.bCL + Tc);


    % Dynamics 
    dy(1)  = -p.dL*L*ePath; %dLdt
    dy(2)  = p.lL*(p.dL)*L*(1+ePath) - iPath ; %dAdt
    dy(3)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;%dThdt
    dy(4)  = lPath + pPath*qPath - tPath; %dTcdt
    dy(5)  = vPath - zPath*uPath; %dTrdt
    dy(6)  = nPath + oPath -p.lC*pPath*qPath -p.lH*sPath*rPath - p.lR*zPath*uPath - wPath; %dIdt

    dcdt = [dy(1), dy(2), dy(3), dy(4), dy(5), dy(6)]';
end
