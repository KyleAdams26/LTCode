function p = parameters()
    % Parameters of Liver Transplant model 
    % We considered parameters and their bounds.
    % The credit of this code is for Dr. Jaimit Parikh.
    % Further credit goes to Dr. Mahya Aghaee for modification of Dr. Jaimit
    % Parikh's code.
    
    % parameters: 
    p.lL = 0.00000452; %1
    p.dA = 0.0833; %3

    p.sR = 0.107; %4
    p.dR = 0.0658; %5
    p.aIR = 0.625; %6
    p.bIR = 0.00833; %7
    
    p.aCI = 0.36; %8
    p.bCI = 352; %9
    p.aHI = 70.7; %10
    p.bHI = 99.7; %11
    p.lC = 0.001; %12
    p.gC = 2.08; %13
    p.KC = 598; %14
    p.aIC = 2;
    p.bIC = .178; %15
    p.lH = 0.0000063; %16
    p.gH = 1.51; %17
    p.KH = 422; %18
    p.aIH = 2;
    p.bIH = .178; %19
    p.lR = 0.000000667; %20
    p.dI = 166; %21
    
    p.aHC = 1; %22
    p.bHC = 35; %23
    p.dC = 0.585; %24
    
    p.aAH = 0.0000261; %25
    p.bAH = 4; %26
    p.aRA = 0.4; %27
    p.bRA = 20; %28
    p.aIRA = 2; %29
    p.bIRA = .356; %30
    p.dH = 0.333; %31

    p.dL = 0.005; %32
    p.aCL = 10; %33
    p.bCL = 200; %34

end
