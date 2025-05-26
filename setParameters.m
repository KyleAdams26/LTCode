function p = setParameters()
    % Parameters of Liver Transplant model (35 parameters)
    % This code is based on code originally written by Jaimit Parikh
    % in 2023 and modified by Mahya Aghaee in 2024
    
    % parameters: 
    p.dL   = 0.005;       % 1
    p.aCL  = 10;          % 2
    p.bCL  = 200;         % 3
    p.lL   = 4.52e-8;     % 4 %may try e-8 in another run
    p.dA   = 0.0833;      % 5
    p.aAH  = 0.0000261;   % 6
    p.bAH  = 4;           % 7
    p.aRA  = 0.4;         % 8
    p.bRA  = 20;          % 9
    p.aIRA = 2;           %10
    p.bIRA = 0.356;       %11
    p.gH   = 1.51;        %12
    p.KH   = 282;         %13
    p.aIH  = 2;           %14
    p.bIH  = 0.178;       %15
    p.dH   = 0.333;       %16
    p.aHC  = 1;           %17
    p.bHC  = 35;          %18
    p.gC   = 2.08;        %19
    p.KC   = 399;         %20
    p.aIC  = 2;           %21
    p.bIC  = 0.178;       %22
    p.dC   = 0.585;       %23
    p.sR   = 0.107;       %24
    p.dR   = 0.0658;      %25
    p.aIR  = 0.625;       %26
    p.bIR  = 0.00833;     %27
    p.aCI  = 0.36;        %28
    p.bCI  = 352;         %29
    p.aHI  = 70.7;        %30
    p.bHI  = 99.7;        %31
    p.lC   = 0.001;       %32
    p.lH   = 0.0000063;   %33
    p.lR   = 0.000000667; %34
    p.dI   = 166;         %35
 
end
