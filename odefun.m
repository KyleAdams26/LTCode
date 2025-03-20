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
    A = c(1); % Alloantigen
    Tr = c(2); % Treg cells
    I = c(3); % IL-2
    Tc = c(4); % Tc cells  
    Th = c(5); % Helper T cells
    L = c(6); % Liver cells

    %     # t : time (int)
    %     # x : 35 parameter values

        %% -- The rate change of the 6 populations (Dynamics of the system) -- %%
        % Paths for the dynamics %
   
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
    dy(1)  = p.lL*(p.dL)*L*(1+ePath) - iPath ; %dAdt
    dy(2)  = vPath - zPath*uPath; %dTrdt
    dy(3)  = nPath + oPath -p.lC*pPath*qPath -p.lH*sPath*rPath - p.lR*zPath*uPath - wPath;
    dy(4)  = lPath + pPath*qPath - tPath; %dTcdt
    dy(5)  = gPath*(1 - jPath*xPath) + sPath*rPath - yPath;
    dy(6) = -p.dL*L*ePath; %dLdt

    dcdt = [dy(1), dy(2), dy(3), dy(4), dy(5), dy(6)]';
end
