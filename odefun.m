function dcdt = odefun(~, c, p) %QC'd
    % odefun - system of equations of the Liver Transplant model
    % Inputs:
    %   ~ - time
    %   c - state variables
    %   p - parameter structure
    % Outputs:
    % dcdt - Derivatives of the state variables
    % The credit of this code is for Dr. Jaimit Parikh 
    % This code was originally modified by Dr. Mahya Aghaee.
    % Kyle Adams modified this code to be suitable for the Liver Transplant
    % Project.
 
    % # c : levels of 6 biological factors at t. 
    L = c(1); % Liver cells
    A = c(2); % APCs
    Th = c(3); % Helper T cells
    Tc = c(4); % Tc cells  
    Tr = c(5); % T reg cells
    I = c(6); % IL-2

    %     # t : time (int)
    %     # x : 35 parameter values

        %% -- The rate change of the 6 populations (Dynamics of the system) -- %%
        % Paths for the dynamics %
   
    ePath = p.aCL*Tc/(p.bCL + Tc);
    gPath = p.aAH*A/(p.bAH + A);
    hPath = p.dL*L;
    iPath = p.dA*A;
    jPath = p.aRA*Tr/(p.bRA + Tr);
    lPath = p.aHC*Th/(p.bHC + Th);
    nPath = p.aCI*Tc/(p.bCI +Tc);
    oPath = p.aHI*Th/(p.bHI + Th);
    pPath = p.gC*Tc*(1 - Tc/p.KC);
    qPath = p.aIC*I/(p.bIC + I);
    rPath = p.aIH*I/(p.bIH + I);
    sPath = p.gH*Th*(1 - Th/p.KH);
    tPath = p.dC*Tc;
    uPath = 1 - p.aIR*I/(p.bIR + I);
    vPath = p.sR;
    wPath = p.dI*I;
    xPath = 1 + p.aIRA*I/(p.bIRA + I);
    yPath = p.dH*Th;
    zPath = p.dR*Tr;



    % Dynamics 
    dy(1) = -hPath*ePath; %dLdt
    dy(2)  = p.lL*(p.dL)*L*(1+ePath) - iPath ; %dAdt
    dy(3)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;%dThdt
    dy(4)  = lPath + pPath*qPath - tPath; %dTcdt
    dy(5)  = vPath - zPath*uPath; %dTrdt
    dy(6)  = nPath + oPath -p.lC*pPath*qPath -p.lH*sPath*rPath - p.lR*zPath*uPath - wPath;

    dcdt = [dy(1), dy(2), dy(3), dy(4), dy(5), dy(6)]';
end