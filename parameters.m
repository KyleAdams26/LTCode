function p = parameters()
    % Parameters of Liver Transplant model 
    % We considered parameters and their bounds.
    % The credit of this code is for Dr. Jaimit Parikh.
    % Further credit goes to Dr. Mahya Aghaee for modification of Dr. Jaimit
    % Parikh's code.
    
    % parameters: 
    p.lL = 0.00000000452; %1
    p.dA = 0.0833; %2

    p.sR = 0.107; %3
    p.dR = 0.0658; %4
    p.aIR = 0.625; %5
    p.bIR = 0.00833; %6
    
    p.aCI = 0.36; %7
    p.bCI = 352; %8
    p.aHI = 70.7; %9
    p.bHI = 99.7; %10
    p.lC = 0.001; %11
    p.gC = 2.08; %12
    p.KC = 527; %13
    p.aIC = 2; %14
    p.bIC = .178; %15
    p.lH = 0.0000063; %16
    p.gH = 1.51; %17
    p.KH = 593; %18
    p.aIH = 2; %19
    p.bIH = .178; %20
    p.lR = 0.000000667; %21
    p.dI = 166; %22
    
    p.aHC = 1; %23
    p.bHC = 35; %24
    p.dC = 0.585; %25
    
    p.aAH = 0.0000261; %26
    p.bAH = 4; %27
    p.aRA = 0.4; %28
    p.bRA = 20; %29
    p.aIRA = 2; %30
    p.bIRA = .356; %31
    p.dH = 0.333; %32

    p.dL = 0.005; %33
    p.aCL = 10; %34
    p.bCL = 200; %35

end