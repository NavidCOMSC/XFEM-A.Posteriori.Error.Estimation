% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [Sxx,Sxy,Syy,Svm] = elemStress(omega,DISPLACEMENT)
% This function calculates the stress distribution within each element from
% the nodal displacements.  The stresses in the xx, yy and xy directions 
% are calculated.

global CHI CONNEC DOMAIN CRACK MAT NODES PLOT PHI PSI XYZ ZETA

nXElem   = DOMAIN(1);                                                       % Number of elements in the x-direction
nYElem   = DOMAIN(2);                                                       % Number of elements in the y-direction
lXElem   = DOMAIN(3);                                                       % Length of elements
Em       = MAT(1);                                                          % Young's modulus for the matrix
vm       = MAT(2);                                                          % Poisson's ratio for the matrix
Ef       = MAT(3);                                                          % Young's modulus for the fiber
vf       = MAT(4);                                                          % Poisson's ratio for the fiber
plane    = MAT(5);                                                          % Plane stress or plane strain
nCT      = size(PHI,2);                                                     % Number of crack tips 
nElem    = nXElem*nYElem;                                                   % Number of elements
nNode    = (nXElem+1)*(nYElem+1);                                           % Number of nodes
stressXX = zeros(nElem,4);                                                  % Create a matrix for storing nodal stress values
stressXY = zeros(nElem,4);                                                  % Create a matrix for storing nodal stress values
stressYY = zeros(nElem,4);                                                  % Create a matrix for storing nodal stress values
stressVM = zeros(nElem,4);                                                  % Create a matrix for storing nodal stress values

m = size(CRACK,1);                                                          % Determine number of data points defining crack
if m > 0
    if nCT == 1
        xCT = CRACK(m,1);                                                   % X-coordinate of crack tip
        yCT = CRACK(m,2);                                                   % Y-coordinate of crack tip
    elseif nCT == 2
        xCT = [CRACK(1,1) CRACK(m,1)];                                      % X-coordinates of crack tips
        yCT = [CRACK(1,2) CRACK(m,2)];                                      % Y-coordinates of crack tips
    end
end

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

%gp = [-1.0 -1.0; 1.0 -1.0; 1.0  1.0;-1.0  1.0];                            % Gauss points defining the nodes of the elements
%gp = [0, 0];
gp = [-0.99 -0.99; 0.99 -0.99; 0.99  0.99;-0.99  0.99];

% Calculate stress at each gauss point
for iElem = 1:(nXElem*nYElem)
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element
    %disp(NN)
    CTN = nnz(NN(:,4));                                                     % Number of nodes with crack tip enrichment    
    HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    IEN = nnz(NN(:,30));                                                    % Number of inclusion nodes    
    NEN = HEN+CTN+IEN;                                                      % Number of crack tip enriched nodes

    % Elemental displacement = [u1;v1;u2;v2;u3;v3;u4;v4]
    U = [DISPLACEMENT(2*N1-1) DISPLACEMENT(2*N1) DISPLACEMENT(2*N2-1) DISPLACEMENT(2*N2)...
         DISPLACEMENT(2*N3-1) DISPLACEMENT(2*N3) DISPLACEMENT(2*N4-1) DISPLACEMENT(2*N4)];

    if (NEN == 0)                                                           % Unenriched element
        for iGP = 1:length(gp)
            
            xi = gp(iGP,1); eta = gp(iGP,2);                                % Current gauss points

            Nx = 2/lXElem*1/4*[-(1-eta);1-eta;1+eta;-(1+eta)];
            Ny = 2/lXElem*1/4*[-(1-xi);-(1+xi);1+xi;1-xi];

            Zgp = ZETA(NN(iGP));                                            % Material level set at current node
            Cgp = CHI(NN(iGP));                                             % Void level set at current node
            
            Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0   Nx(4)   0;...
                    0   Ny(1)   0   Ny(2)   0   Ny(3)   0   Ny(4);...
                  Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];
            
            if Cgp >= 0
                if Zgp == 0, Zgp = setdiff(ZETA(NN(:,1)),0); end
                if Zgp > 0, C = Cm; else C = Cf; end
                stress = C*Bu*U';
            else
                stress = zeros(1,3);
            end
            
            stressXX(iElem,iGP) = stress(1);
            stressYY(iElem,iGP) = stress(2);
            stressXY(iElem,iGP) = stress(3);
            stressVM(iElem,iGP) = sqrt(stress(1)^2+stress(2)^2-stress(1)*stress(2)+3*stress(3)^2);            
        end
    else
        Uenr = [];
        for iGP = 1:length(gp)
            xi = gp(iGP,1); eta = gp(iGP,2);

            N  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
            Nx = 2/lXElem*1/4*[-(1-eta);1-eta;1+eta;-(1+eta)];
            Ny = 2/lXElem*1/4*[-(1-xi);-(1+xi);1+xi;1-xi];

            Xgp = XYZ(NN(iGP,1),2);                                         % The global X for the current gauss point
            Ygp = XYZ(NN(iGP,1),3);                                         % The global Y for the current gauss point
            Zgp = ZETA(NN(iGP));                                            % Material level set at current gauss point                     
            
            Benr = [];
            Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0   Nx(4)   0;...
                    0   Ny(1)   0   Ny(2)   0   Ny(3)   0   Ny(4);...
                  Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];

            iB = 1; iLoc = 1;  
            for iN = 1:4
                if NN(iN,2) ~= 0
                    psi1 = PSI(N1);                                         % Psi level set value at node 1
                    psi2 = PSI(N2);                                         % Psi level set value at node 2
                    psi3 = PSI(N3);                                         % Psi level set value at node 3
                    psi4 = PSI(N4);                                         % Psi level set value at node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
                    if psi == 1e-6, psi = 0; end
                    Hgp  = sign(psi);
                    
                    if Hgp == 0, Hgp = sign(nonzeros([psi1 psi2 psi3 psi4])); end                    
                    if length(Hgp) > 1, Hgp = sign(Hgp(1)); end
                    Hi = NN(iN,3);
                    H  = Hgp-Hi;
                    
                    Ba = [Nx(iN)*H     0;
                              0    Ny(iN)*H;
                          Ny(iN)*H Nx(iN)*H];
                    Benr(:,iB:(iB+1)) = Ba;
                    iB = iB+2;

                    if iGP == 1
                        Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NN(iN,2)-1) DISPLACEMENT(2*NN(iN,2))];
                        iLoc = iLoc+2;
                    end
                elseif NN(iN,4) ~= 0
                    if nCT == 1
                        X     = Xgp-xCT;                                    % Horizontal distance from crack tip to gauss point
                        Y     = Ygp-yCT;                                    % Vertical distance from crack tip to gauss point
                        CCS   = [cos(omega) sin(omega);-sin(omega) cos(omega)];
                        XYloc = CCS*[X Y]';                                 % Change to crack tip coordinates
                        r     = sqrt(XYloc(1)^2+XYloc(2)^2);                % Radius from crack tip to current gauss point
                        if r < 0.001*lXElem; r = 0.05*lXElem; end
                        theta = atan2(XYloc(2),XYloc(1));                   % Angle from crack tip to current gauss point
                    elseif nCT == 2
                        X1  = Xgp-xCT(1);
                        Y1  = Ygp-yCT(1);
                        X2  = Xgp-xCT(2);
                        Y2  = Ygp-yCT(2);
                        CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
                        XY1 = CCS*[X1 Y1]';
                        CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];
                        XY2 = CCS*[X2 Y2]';
                        r1  = sqrt(XY1(1)^2+XY1(2)^2);                      % Radius from crack tip to current gauss point
                        r2  = sqrt(XY2(1)^2+XY2(2)^2);
                        if r1 > r2
                            r = r2; theta = atan2(XY2(2),XY2(1));
                            CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];                            
                        elseif r2 > r1
                            r = r1; theta = atan2(XY1(2),XY1(1));
                            CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
                        end
                        if r < 0.001*lXElem; r = 0.05*lXElem; end     
                    end

                    c = 1/2/sqrt(r); ct = CCS(1,1); st = CCS(1,2);          % Constants
                    
                    if NN(iN,12) == 0                                       % Crack tip enrichment
                        a1gp  = sqrt(r)*sin(theta/2);                       % Node 1 crack tip enrichment value
                        a2gp  = sqrt(r)*cos(theta/2);                       % Node 2 crack tip enrichment value
                        a3gp  = sqrt(r)*sin(theta)*sin(theta/2);            % Node 3 crack tip enrichment value
                        a4gp  = sqrt(r)*sin(theta)*cos(theta/2);            % Node 4 crack tip enrichment value

                        a1  = a1gp-NN(iN,5); a2  = a2gp-NN(iN,7);           % Shifted crack tip enrichment values
                        a3  = a3gp-NN(iN,9); a4  = a4gp-NN(iN,11);          % Shifted crack tip enrichment values

                        % Derivative of crack tip enrichment functions with respect to X
                        Px = c*[-sin(theta/2)*ct              + cos(theta/2)*-st;...
                                 cos(theta/2)*ct              + sin(theta/2)*-st;...
                                -sin(3*theta/2)*sin(theta)*ct + (sin(theta/2)+sin(3*theta/2)*cos(theta))*-st;...
                                -cos(3*theta/2)*sin(theta)*ct + (cos(theta/2)+cos(3*theta/2)*cos(theta))*-st];

                        % Derivative of crack tip enrichment functions with respect to Y
                        Py = c*[-sin(theta/2)*st              + cos(theta/2)*ct;...
                                 cos(theta/2)*st              + sin(theta/2)*ct;...
                                -sin(3*theta/2)*sin(theta)*st + (sin(theta/2)+sin(3*theta/2)*cos(theta))*ct;...
                                -cos(3*theta/2)*sin(theta)*st + (cos(theta/2)+cos(3*theta/2)*cos(theta))*ct];

                        B1 = [Nx(iN)*a1+N(iN)*Px(1)          0;
                                       0            Ny(iN)*a1+N(iN)*Py(1);
                              Ny(iN)*a1+N(iN)*Py(1) Nx(iN)*a1+N(iN)*Px(1)];

                        B2 = [Nx(iN)*a2+N(iN)*Px(2)          0;
                                     0            Ny(iN)*a2+N(iN)*Py(2);
                              Ny(iN)*a2+N(iN)*Py(2) Nx(iN)*a2+N(iN)*Px(2)];

                        B3 = [Nx(iN)*a3+N(iN)*Px(3)          0;
                                     0            Ny(iN)*a3+N(iN)*Py(3);
                              Ny(iN)*a3+N(iN)*Py(3) Nx(iN)*a3+N(iN)*Px(3)];

                        B4 = [Nx(iN)*a4+N(iN)*Px(4)          0;
                                     0            Ny(iN)*a4+N(iN)*Py(4);
                              Ny(iN)*a4+N(iN)*Py(4) Nx(iN)*a4+N(iN)*Px(4)];

                        Bb = [B1 B2 B3 B4];
                        Benr(:,iB:(iB+7)) = Bb;
                        iB = iB+8;

                        if iGP == 1
                            index = 3;
                            for iAlpha = 1:4
                                Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NN(iN,iAlpha+index)-1) DISPLACEMENT(2*NN(iN,iAlpha+index))];
                                index = index+1;
                                iLoc  = iLoc+2;
                            end
                        end
                    else                                                    % Bimaterial crack enrichment
                        G1 = Em/2/(1+vm);                                   % Shear modulus for the matrix
                        G2 = Ef/2/(1+vf);                                   % Shear modulus for the fiber                        
                        
                        % Kosolov constant
                        if plane == 1                                       % Plane stress
                            k1 = (3-vm)/(1+vm);
                            k2 = (3-vf)/(1+vf);
                        elseif plane == 2                                   % Plane strain
                            k1 = 3-4*vm;
                            k2 = 3-4*vf;
                        end
                        
                        b = (G1*(k2-1)-G2*(k1-1))/(G1*(k2+1)+G2*(k1+1));    % Second Dundur's parameter
                        e = 1/(2*pi)*log((1-b)/(1+b));                      % Material constant 
                        
                        % Common variables
                        sr  = sqrt(r);
                        st  = sin(theta);
                        st2 = sin(theta/2);
                        ct2 = cos(theta/2);
                       
                        a1gp  = sr*cos(e*log(r))*exp(-e*theta)*st2;         % Alpha 1 crack tip enrichment value
                        a2gp  = sr*cos(e*log(r))*exp(-e*theta)*ct2;         % Alpha 2 crack tip enrichment value
                        a3gp  = sr*cos(e*log(r))*exp(e*theta)*st2;          % Alpha 3 crack tip enrichment value
                        a4gp  = sr*cos(e*log(r))*exp(e*theta)*ct2;          % Alpha 4 crack tip enrichment value
                        a5gp  = sr*cos(e*log(r))*exp(e*theta)*st2*st;       % Alpha 5 crack tip enrichment value
                        a6gp  = sr*cos(e*log(r))*exp(e*theta)*ct2*st;       % Alpha 6 crack tip enrichment value
                        a7gp  = sr*sin(e*log(r))*exp(-e*theta)*st2;         % Alpha 7 crack tip enrichment value
                        a8gp  = sr*sin(e*log(r))*exp(-e*theta)*ct2;         % Alpha 8 crack tip enrichment value
                        a9gp  = sr*sin(e*log(r))*exp(e*theta)*st2;          % Alpha 9 crack tip enrichment value
                        a10gp = sr*sin(e*log(r))*exp(e*theta)*ct2;          % Alpha 10 crack tip enrichment value
                        a11gp = sr*sin(e*log(r))*exp(e*theta)*st2*st;       % Alpha 11 crack tip enrichment value
                        a12gp = sr*sin(e*log(r))*exp(e*theta)*ct2*st;       % Alpha 12 crack tip enrichment value
                        
                        a = [a1gp-NN(iN,5);...                              % Shifted alpha 1 enrichment value
                             a2gp-NN(iN,7);...                              % Shifted alpha 2 enrichment value
                             a3gp-NN(iN,9);...                              % Shifted alpha 3 enrichment value
                             a4gp-NN(iN,11);...                             % Shifted alpha 4 enrichment value
                             a5gp-NN(iN,13);...                             % Shifted alpha 5 enrichment value
                             a6gp-NN(iN,15);...                             % Shifted alpha 6 enrichment value
                             a7gp-NN(iN,17);...                             % Shifted alpha 7 enrichment value
                             a8gp-NN(iN,19);...                             % Shifted alpha 8 enrichment value
                             a9gp-NN(iN,21);...                             % Shifted alpha 9 enrichment value
                             a10gp-NN(iN,23);...                            % Shifted alpha 10 enrichment value
                             a11gp-NN(iN,25);...                            % Shifted alpha 11 enrichment value
                             a12gp-NN(iN,27)];                              % Shifted alpha 12 enrichment value
                        
                        % Derivative of bimaterial crack tip enrichment functions with respect to x1 (crack tip coordinate system)
                        px = c*[-exp(-e*theta)*sin(theta/2)*(cos(e*log(r))+2*e*sin(e*log(r)-theta));...
                                 exp(-e*theta)*cos(theta/2)*(cos(e*log(r))-2*e*sin(e*log(r)-theta));...
                                -exp(e*theta)*sin(theta/2)*(cos(e*log(r))+2*e*sin(e*log(r)+theta));...
                                 exp(e*theta)*cos(theta/2)*(cos(e*log(r))-2*e*sin(e*log(r)+theta));...
                                -exp(e*theta)*sin(theta)*(cos(e*log(r))*sin(3*theta/2)+2*e*sin(e*log(r)+theta)*sin(theta/2));...
                                -exp(e*theta)*sin(theta)*(cos(e*log(r))*cos(3*theta/2)+2*e*sin(e*log(r)+theta)*cos(theta/2));...
                                 exp(-e*theta)*sin(theta/2)*(-sin(e*log(r))+2*e*cos(e*log(r)-theta));...
                                 exp(-e*theta)*cos(theta/2)*(sin(e*log(r))+2*e*cos(e*log(r)-theta));...
                                 exp(e*theta)*sin(theta/2)*(-sin(e*log(r))+2*e*cos(e*log(r)+theta));...
                                 exp(e*theta)*cos(theta/2)*(sin(e*log(r))+2*e*cos(e*log(r)+theta));...
                                 exp(e*theta)*sin(theta)*(-sin(e*log(r))*sin(3*theta/2)+2*e*cos(e*log(r)+theta)*sin(theta/2));...
                                 exp(e*theta)*sin(theta)*(-sin(e*log(r))*cos(3*theta/2)+2*e*cos(e*log(r)+theta)*cos(theta/2))];

                        % Derivative of bimaterial crack tip enrichment functions with respect to x2 (crack tip coordinate system)
                        py = c*[exp(-e*theta)*(cos(e*log(r))*cos(theta/2)-2*e*cos(e*log(r)-theta)*sin(theta/2));...
                                exp(-e*theta)*(cos(e*log(r))*sin(theta/2)-2*e*cos(e*log(r)-theta)*cos(theta/2));...
                                exp(e*theta)*(cos(e*log(r))*cos(theta/2)+2*e*cos(e*log(r)+theta)*sin(theta/2));...
                                exp(e*theta)*(cos(e*log(r))*sin(theta/2)+2*e*cos(e*log(r)+theta)*cos(theta/2));...
                                exp(e*theta)*(cos(e*log(r))*(sin(theta/2)+sin(3*theta/2)*cos(theta))+2*e*cos(e*log(r)+theta)*sin(theta/2)*sin(theta));...
                                exp(e*theta)*(cos(e*log(r))*(cos(theta/2)+cos(3*theta/2)*cos(theta))+2*e*cos(e*log(r)+theta)*cos(theta/2)*sin(theta));...
                                exp(-e*theta)*(sin(e*log(r))*cos(theta/2)-2*e*sin(e*log(r)-theta)*sin(theta/2));...
                                exp(-e*theta)*(sin(e*log(r))*sin(theta/2)-2*e*sin(e*log(r)-theta)*cos(theta/2));...
                                exp(e*theta)*(sin(e*log(r))*cos(theta/2)+2*e*sin(e*log(r)+theta)*sin(theta/2));...
                                exp(e*theta)*(sin(e*log(r))*sin(theta/2)+2*e*sin(e*log(r)+theta)*cos(theta/2));...
                                exp(e*theta)*(sin(e*log(r))*(sin(theta/2)+sin(3*theta/2)*cos(theta))+2*e*sin(e*log(r)+theta)*sin(theta/2)*sin(theta));...
                                exp(e*theta)*(sin(e*log(r))*(cos(theta/2)+cos(3*theta/2)*cos(theta))+2*e*sin(e*log(r)+theta)*cos(theta/2)*sin(theta))];
                        
                        % Derivative of bimaterial crack tip enrichment functions with respect to X
                        Px = [px(1)*ct+py(1)*-st;...
                              px(2)*ct+py(2)*-st;...
                              px(3)*ct+py(3)*-st;...
                              px(4)*ct+py(4)*-st;...
                              px(5)*ct+py(5)*-st;...
                              px(6)*ct+py(6)*-st;...
                              px(7)*ct+py(7)*-st;...
                              px(8)*ct+py(8)*-st;...
                              px(9)*ct+py(9)*-st;...
                              px(10)*ct+py(10)*-st;...
                              px(11)*ct+py(11)*-st;...
                              px(12)*ct+py(12)*-st];

                        % Derivative of bimaterial crack tip enrichment functions with respect to Y
                        Py = [px(1)*st+py(1)*ct;...
                              px(2)*st+py(2)*ct;...
                              px(3)*st+py(3)*ct;...
                              px(4)*st+py(4)*ct;...
                              px(5)*st+py(5)*ct;...
                              px(6)*st+py(6)*ct;...
                              px(7)*st+py(7)*ct;...
                              px(8)*st+py(8)*ct;...
                              px(9)*st+py(9)*ct;...
                              px(10)*st+py(10)*ct;...
                              px(11)*st+py(11)*ct;...
                              px(12)*st+py(12)*ct];

                        Bb = zeros(3,24);
                        for iA = 1:12
                            Balpha = [Nx(iN)*a(iA)+N(iN)*Px(iA)             0;...
                                                  0             Nx(iN)*a(iA)+N(iN)*Py(iA);...
                                      Nx(iN)*a(iA)+N(iN)*Py(iA) Nx(iN)*a(iA)+N(iN)*Px(iA)];
                            Bb(:,(2*iA-1):2*iA) = Balpha;
                        end
                        Benr(:,iB:(iB+23)) = Bb;
                        iB = iB+24;
                        
                        if iGP == 1
                            index = 3;
                            for iAlpha = 1:12
                                Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NN(iN,iAlpha+index)-1) DISPLACEMENT(2*NN(iN,iAlpha+index))];
                                index = index+1;
                                iLoc  = iLoc+2;
                            end
                        end                        
    
                    end
                end
                
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
            stress = C*[Bu Benr]*[U Uenr]';
            stressXX(iElem,iGP) = stress(1);
            stressYY(iElem,iGP) = stress(2);
            stressXY(iElem,iGP) = stress(3);
            stressVM(iElem,iGP) = sqrt(stress(1)^2+stress(2)^2-stress(1)*stress(2)+3*stress(3)^2);
        end
    end
end

% Average nodal stress values
if PLOT(4,2) == 1                                                           % Average nodal stress values
    Sxx = zeros(nNode,2); Syy = Sxx; Sxy = Sxx; Svm = Sxx;

    % Construct stress vectors
    for iE = 1:(nXElem*nYElem)
        for iN = 1:4
            nNode = CONNEC(iE,iN+1);
            Sxx(nNode,:) = Sxx(nNode,:) + [stressXX(iE,iN) 1];
            Sxy(nNode,:) = Sxy(nNode,:) + [stressXY(iE,iN) 1];
            Syy(nNode,:) = Syy(nNode,:) + [stressYY(iE,iN) 1];
            Svm(nNode,:) = Svm(nNode,:) + [stressVM(iE,iN) 1];
        end
    end
    
    % Average nodal stress values
    Sxx(:,1) = Sxx(:,1)./Sxx(:,2); Sxx(:,2) = [];
    Sxy(:,1) = Sxy(:,1)./Sxy(:,2); Sxy(:,2) = [];
    Syy(:,1) = Syy(:,1)./Syy(:,2); Syy(:,2) = [];
    Svm(:,1) = Svm(:,1)./Svm(:,2); Svm(:,2) = [];
else                                                                        % Do not average nodal stress values
    Sxx = stressXX;
    Sxy = stressXY;
    Syy = stressYY;
    Svm = stressVM;
end