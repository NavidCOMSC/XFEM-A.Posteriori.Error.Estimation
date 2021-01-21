% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

clear all; close all; format compact; tic; global CHI CRACK NODES ZETA CONNEC

%%%% Pre-Processing %%%%%
inputQuasiStatic;                                                           % Define the geometry, materials, discontinuities

iter = numIterations;

for i = 1:iter
    
    %%%%% Processing %%%%%
    if i == 1, connectivity; pHDOF = []; else NODES(:,4:29) = 0; end        % Define connectivity
    omega           = levelSet(i);                                          % Create phi and psi, define enriched elements
    [DOF,DISP]      = calcDOF;                                              % Total degrees of freedom
    [updElem,IElem] = enrElem(i,pHDOF);                                     % Find enriched elements, inclusion elements
    
    if i == 1,[globalK, globalF] = stiffnessMatrix(omega,DOF,iter,updElem,IElem);    % Construct global stiffness matrix        
    else globalK = updateStiffness(globalK,omega,DOF,updElem,pHDOF); end    % Update the global stiffness matrix
    
    %disp(globalF)
    globalF         = forceVector(DOF,i,globalF);                           % Construct global force vector
    %disp(globalF')
    freeDOF         = boundaryCond(DOF);                                    % Solve for the degrees of freedom
    
    %globalF
    
    %DISP(freeDOF,:) = globalK(freeDOF,freeDOF)\globalF(freeDOF,:);         % Find the nodal displacement values
    temp=sqrt(size(NODES,1));
    %freeDOF=2*temp+1:2*temp*temp;
    DISP(freeDOF) = globalK(freeDOF,freeDOF)\globalF(freeDOF);             % Find the nodal displacement values
    pause
    %DISP(freeDOF) = pinv(full(globalK(freeDOF,freeDOF)))*globalF(freeDOF);    
    %%% Post-Processing %%%%%
    if i == iter, plotMain(omega,DISP); end                                 % Make plots
    if isempty(CRACK) == 0
        pHDOF    = 2*max(NODES(:,2));                                       % Maximum constant DOF at current iteration 
        [KI,KII] = JIntegral(omega,DISP);                                   % Calculate the stress intensity factors
        exit     = growCrack(KI,KII,omega);                                 % Advance crack for quasi-static growth
       
        if strcmp(exit,'YES') == 1
            disp('WARNING: No crack growth, iterations exited early.')
            plotMain(omega,DISP); break 
        end
    end
    %[uxL2Norm,uyL2Norm,L2Norm] =GetL2Norm(DISP);
    %sqrt(uxL2Norm)
    %sqrt(uyL2Norm)
    %sqrt(L2Norm)
    
    [EnergyNorm, EnergyNormElem] = GetEnergyNorm(DISP);
    sqrt(EnergyNorm)
    
     figure
     y=linspace(0,1,1001);
     numinterval = 11;
     [~,uy] = getExactSol(0.5,y,MAT,FORCE,INC);
     uyap = DISP(2:2*(numinterval+1):(numinterval+1)^2*2);
     plot(y,uy,'-b',0:1/numinterval:1,uyap,'-g')
     
    [DUsU, SSU ] = GetSueprconvergentUpdU(DISP );
    [DUsD, SSD ] = GetSueprconvergentDowndU(DISP );
    
    %DUsU
    %DUsD
    
    
    
    %SpConTotal = [DUsU, DUsD, SSU, SSD];
    %SS
    
    [AvEU, AvSU ] = avmatrix(DUsU, SSU );
    [AvED, AvSD ] = avmatrix(DUsD, SSD );
    
    
    %AvEU
    %AvED
    
    AvEU((numinterval+1)*(numinterval-1)/2+1:(numinterval+1)*(numinterval+1)/2, :) = ...
        AvED((numinterval+1)*(numinterval-1)/2+1:(numinterval+1)*(numinterval+1)/2, :);

    AvE = AvEU;
    AvS = AvSU;
    
    DUs = (DUsU + DUsD)./2;
    SS = (SSU + SSD)./2;
    
    %disp(AvE)
    %size(AvS)
    %AvS
    Posteriori = Getaposteriori( AvS, DISP );
    Posteriori
    Posteriori./EnergyNormElem
    
    disp(['Iteration ',num2str(i),' completed. Elapsed time is ',num2str(toc,'%0.4f'),'.'])    
end