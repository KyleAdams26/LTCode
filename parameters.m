function p = parameters()
    % Parameters of Liver Transplant model 
    % We considered parameters and their bounds.
    % The credit of this code is for Dr. Jaimit Parikh.
    % Further credit goes to Dr. Mahya Aghaee for modification of Dr. Jaimit
    % Parikh's code.
    
    % parameters: 
    p.aLA = 1; %1
    p.bLA = 5; %2
    p.dA = 0.083333; %3

    p.sR = 0.05; %4
    p.dR = 0.0658; %5
    p.aIR = 0.608; %6
    p.bIR = 0.00833333; %7
    
    p.aCI = 0.36; %8
    p.bCI = 350; %9
    p.aHI = 70.7; %10
    p.bHI = 100; %11
    p.lC = 0.005; %12
    p.gC = 2.079; %13
    p.KC = 598.2; %14
    p.bIC = .178; %15
    p.lH = 0.000005999988; %16
    p.gH = 1.512; %17
    p.KH = 422.4; %18
    p.bIH = 1.78; %19
    p.lR = 0.000003333332; %20
    p.dI = 166.355; %21
    
    p.aHC = 1; %22
    p.bHC = 352; %23
    p.dC = 0.5853658537; %24
    
    p.aAH = 0.0000261; %25
    p.bAH = 4; %26
    p.aRA = 0.4; %27
    p.bRA = 20; %28
    p.aIRA = 2; %29
    p.bIRA = 356; %30
    p.dH = 0.3333; %31

    p.dL = 1/200; %32
    p.aCL = 10; %33
    p.bCL = 200; %34

end