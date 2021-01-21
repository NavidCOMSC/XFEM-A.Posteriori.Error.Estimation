function [DUs, SS ] = GetSueprconvergentUpdU(DISPLACEMENT )

global CONNEC DOMAIN MAT NODES ZETA

nXElem   = DOMAIN(1);                                                       % Number of elements in the x-direction
nYElem   = DOMAIN(2);                                                       % Number of elements in the y-direction
lXElem   = DOMAIN(3);                                                       % Length of elements
Em       = MAT(1);                                                          % Young's modulus for the matrix
vm       = MAT(2);                                                          % Poisson's ratio for the matrix
Ef       = MAT(3);                                                          % Young's modulus for the fiber
vf       = MAT(4);                                                          % Poisson's ratio for the fiber
plane    = MAT(5);                                                          % Plane stress or plane strain
%nCT      = size(PHI,2);                                                     % Number of crack tips 
nElem    = nXElem*nYElem;                                                   % Number of elements
%nNode    = (nXElem+1)*(nYElem+1);                                           % Number of nodes
DUs = zeros(nElem,3);                                                      % Create a matrix for storing nodal strain values in x coordinate
SS = zeros(nElem,3);
%DUyy = zeros(nElem,3);                                                      % Create a matrix for storing nodal strain values in y coordinate
%DUxy = zeros(nElem,3);                                                      % Create a matrix for storing nodal strain values in xy coordinate




% Create elastic constant matrix
if plane == 1                                                               % Plane stress
    C1 = Em/(1-vm^2);                                                       % Constant for elastic constant matrix
    C2 = Em*vm/(1-vm^2);                                                    % Constant for elastic constant matrix
    C3 = Em/2/(1+vm);                                                       % Constant for elastic constant matrix
    Cm = [C1 C2  0;...
          C2 C1  0;...
           0  0 C3];
    C1 = Ef/(1-vf^2);                                                       % Constant for elastic constant matrix
    C2 = Ef*vf/(1-vf^2);                                                    % Constant for elastic constant matrix
    C3 = Ef/2/(1+vf);                                                       % Constant for elastic constant matrix
    Cf = [C1 C2  0;...
          C2 C1  0;...
           0  0 C3];        
elseif plane == 2                                                           % Plane strain
    C1 = Em*(1-vm)/(1+vm)/(1-2*vm);                                         % Constant for elastic constant matrix
    C2 = Em*vm/(1+vm)/(1-2*vm);                                             % Constant for elastic constant matrix
    C3 = Em/2/(1+vm);                                                       % Constant for elastic constant matrix
    Cm = [C1 C2  0;...
          C2 C1  0;...
           0  0 C3];
    C1 = Ef*(1-vf)/(1+vf)/(1-2*vf);                                         % Constant for elastic constant matrix
    C2 = Ef*vf/(1+vf)/(1-2*vf);                                             % Constant for elastic constant matrix
    C3 = Ef/2/(1+vf);                                                       % Constant for elastic constant matrix
    Cf = [C1 C2  0;...
          C2 C1  0;...
           0  0 C3];    
end


%gp = [-1 -1; 1 -1; 1  1;-1  1];                                             %Gauss points defining the nodes of elements
gp = [0,0.000001];

%Calculate strain at each superconvergent points
for iElem = 1:(nXElem*nYElem)
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element
    CTN = nnz(NN(:,4));                                                     % Number of nodes with crack tip enrichment    
    HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    IEN = nnz(NN(:,30));                                                    % Number of inclusion nodes    
    NEN = HEN+CTN+IEN;                                                      % Number of crack tip enriched nodes

    % Elemental displacement = [u1;v1;u2;v2;u3;v3;u4;v4]
    U = [DISPLACEMENT(2*N1-1) DISPLACEMENT(2*N1) DISPLACEMENT(2*N2-1) DISPLACEMENT(2*N2)...
         DISPLACEMENT(2*N3-1) DISPLACEMENT(2*N3) DISPLACEMENT(2*N4-1) DISPLACEMENT(2*N4)];

     if (NEN == 0)                                                           % Unenriched element
        for iGP = 1:size(gp,1)
            
            ZN = [ZETA(N1) ZETA(N2) ZETA(N3) ZETA(N4)];
            xi = gp(iGP,1); eta = gp(iGP,2);                                % Current gauss points
            N = 1/4*[(1-xi).*(1-eta) (1+xi).*(1-eta) (1+xi).*(1+eta) (1-xi).*(1+eta)];

            Nx = 2/lXElem*1/4*[-(1-eta);1-eta;1+eta;-(1+eta)];
            Ny = 2/lXElem*1/4*[-(1-xi);-(1+xi);1+xi;1-xi];

            %Zgp = ZETA(NN(iGP));                                            % Material level set at current node
            Zgp = N(1)*ZN(1)+N(2)*ZN(2)+N(3)*ZN(3)+N(4)*ZN(4);  % Material level set at current gauss point
           
            
            Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0   Nx(4)   0;...
                    0   Ny(1)   0   Ny(2)   0   Ny(3)   0   Ny(4);...
                  Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];
              
             if Zgp == 0, Zgp = setdiff(ZETA(NN(:,1)),0); end
             
             if Zgp > 0, C = Cm; else C = Cf; end
             
             DUs(iElem,:) = U*Bu';
             SS(iElem,:) = U*Bu'*C;
             %DUxx(iElem,iGP) = DU(1);
             %DUyy(iElem,iGP) = DU(2);
             %DUxy(iElem,iGP) = DU(3);
             
             
             
        end
     else
         Uenr = [];
        for iGP = 1:size(gp,1)
            xi = gp(iGP,1); eta = gp(iGP,2);
            ZN = [ZETA(N1) ZETA(N2) ZETA(N3) ZETA(N4)];

            N  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
            Nx = 2/lXElem*1/4*[-(1-eta);1-eta;1+eta;-(1+eta)];
            Ny = 2/lXElem*1/4*[-(1-xi);-(1+xi);1+xi;1-xi];

            %Zgp = ZETA(NN(iGP));
            Zgp = N(1)*ZN(1)+N(2)*ZN(2)+N(3)*ZN(3)+N(4)*ZN(4);  % Material level set at current gauss point
            
            Benr = [];
            Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0   Nx(4)   0;...
                    0   Ny(1)   0   Ny(2)   0   Ny(3)   0   Ny(4);...
                  Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];

            iB = 1; iLoc = 1;  
            for iN = 1:4
                % Inclusion enrichment
                if NN(iN,30) ~= 0
                    zeta = ZETA([N1 N2 N3 N4]);
                    Zm   = dot(N,zeta);
                    Za   = dot(N,abs(zeta));

                    E  = Za-abs(Zm);
                    if E ~= 0
                        Ex = dot(Nx,abs(zeta))-(Zm)/abs(Zm)*dot(Nx,zeta);
                        Ey = dot(Ny,abs(zeta))-(Zm)/abs(Zm)*dot(Ny,zeta);
                        
                        Ba = [Nx(iN)*E+N(iN)*Ex         0;
                                      0         Ny(iN)*E+N(iN)*Ey;
                              Ny(iN)*E+N(iN)*Ey Nx(iN)*E+N(iN)*Ex];
                        Benr(:,iB:(iB+1)) = Ba;
                        iB = iB+2;
                        
                        if iGP == 1
                            Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NN(iN,30)-1) DISPLACEMENT(2*NN(iN,30))];
                            iLoc = iLoc+2;
                        end
                    end
                end               
            end
            if Zgp > 0, C = Cm; else C = Cf; end
        DUs(iElem,:) = [U Uenr]*[Bu Benr]';
        SS(iElem,:) = [U Uenr]*[Bu Benr]'*C;
        %DUxx(iElem,iGP) = DU(1);
        %DUyy(iElem,iGP) = DU(2);
        %DUxy(iElem,iGP) = DU(3);
        
        end
      end
end

