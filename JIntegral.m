% Original version of J-Integral code written by Nguyen Vinh Phu (2006)
% Modified By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [KI KII] = JIntegral(omega,DISPLACEMENT)
% This function calculates the mixed-mode stress intensity factors for a
% cracked body using the domain form of the interaction integrals based on
% the original formulation by Rice.

global CHI CONNEC CRACK DOMAIN MAT NODES PSI XYZ ZETA

nElem  = DOMAIN(1)*DOMAIN(2);                                               % Total number of elements
nNode  = (DOMAIN(1)+1)*(DOMAIN(2)+1);                                       % Total number of nodes
lXElem = DOMAIN(3);                                                         % Length of elements in the x-direction
Em     = MAT(1);                                                            % Young's modulus for the matrix
vm     = MAT(2);                                                            % Poisson's ratio for the matrix
Ef     = MAT(3);                                                            % Young's modulus for the fiber
vf     = MAT(4);                                                            % Poisson's ratio for the fiber
plane  = MAT(5);                                                            % Plane stress or plane strain
Gm     = Em/2/(1+vm);                                                       % Shear modulus for the matrix
Gf     = Ef/2/(1+vf);                                                       % Shear modulus for the fiber

% Create elastic constant matrix
if plane == 1                                                               % Plane stress
    t = MAT(6);                                                             % Plane stress thickness
    
    C1 = Em/(1-2*vm);                                                       % Constant for elastic constant matrix
    C2 = Em*vm/(1-2*vm);                                                    % Constant for elastic constant matrix
    C3 = Em/2/(1+vm);                                                       % Constant for elastic constant matrix
    Cm = t*[C1 C2  0;...
            C2 C1  0;...
             0  0 C3];
    C1 = Ef/(1-2*vf);                                                       % Constant for elastic constant matrix
    C2 = Ef*vf/(1-2*vf);                                                    % Constant for elastic constant matrix
    C3 = Ef/2/(1+vf);                                                       % Constant for elastic constant matrix
    Cf = t*[C1 C2  0;...
            C2 C1  0;...
             0  0 C3];
    km = (3-vm)/(1+vm);
    kf = (3-vf)/(1+vf);
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
    km = 3-4*vm;
    kf = 3-4*vf;
end

nCT = length(omega);                                                        % Number of crack tips in domain

for iJ = 1:nCT

    if iJ == 1
        nPt = size(CRACK,1);                                                % Determine number of data points defining crack
        xCT = CRACK(nPt,1);                                                 % X-coordinate of crack tip
        yCT = CRACK(nPt,2);                                                 % Y-coordinate of crack tip
    else
        xCT = CRACK(1,1);                                                   % X-coordinate of crack tip
        yCT = CRACK(1,2);                                                   % Y-coordinate of crack tip
    end
    
    CCS = [cos(omega(iJ)) sin(omega(iJ));-sin(omega(iJ)) cos(omega(iJ))];
    
    % Geometric predicates for crack geometry
    area  = lXElem*lXElem;                                                  % Elemental area
    
    % Define elements to be used in J-integral
    c      = 3;                                                             % Magnification factor
    radius = c*sqrt(area);                                                  % Search radius
    
    dist = zeros(1,nNode);
    for iN = 1:nNode
        Xn       = XYZ(iN,2);                                               % X-coordinate for the current node
        Yn       = XYZ(iN,3);                                               % Y-coordinate for the current node
        X        = Xn-xCT;                                                  % Horizontal distance from crack tip to current node
        Y        = Yn-yCT;                                                  % Vertical distance from crack tip to current node
        XYloc    = CCS*[X Y]';                                              % Change to crack tip coordinates
        r        = sqrt(XYloc(1)^2+XYloc(2)^2);                             % Radius from crack tip to current gauss point
        dist(iN) = r;                                                       % Store radius value
    end
    
    % Determine elements in J-integral and assign nodal q values
    temp    = dist-radius;                                                  % Determine whether or not the node is outside the search radius
    temp    = temp(CONNEC(:,2:5))';                                         % Build elemental distance vector
    Jdomain = NaN(1,nElem);                                                 % Initialize Jdomain to NaN
    index   = 1;                                                            % Index to track Jdomain
    for i = 1:nElem
        if (min(temp(:,i)) < 0)
            Jdomain(index) = i;
            index = index+1;
        end
    end
    
    Jdomain(isnan(Jdomain)) = [];                                           % Remove the unused locations in Jdomain
    temp  = dist-radius;                                                    % Determine whether or not the node is outside the search radius
    temp  = temp(CONNEC(Jdomain,2:5))';                                     % Build elemental distance vector for elements in J-integral domain
    temp  = (temp<=0);                                                      % Calculate nodes inside/outside search radius
    qNode = temp';                                                          % Store the nodal q values used in J-integral calculations
    
    I = zeros(2,1);                                                         % Interaction integral for combined states 1 and 2
    for iElem = 1:length(Jdomain)
        nElem = Jdomain(iElem);                                             % The global number for the current element
        N1  = CONNEC(nElem,2);                                              % Node 1 for current element
        N2  = CONNEC(nElem,3);                                              % Node 2 for current element
        N3  = CONNEC(nElem,4);                                              % Node 3 for current element
        N4  = CONNEC(nElem,5);                                              % Node 4 for current element
        NN  = NODES([N1 N2 N3 N4]',:);                                      % Nodal data for current element
        CN  = [CHI(N1) CHI(N2) CHI(N3) CHI(N4)];                            % Nodal chi level set values
        CTN = nnz(NN(:,4));                                                 % Number of nodes with crack tip enrichment
        HEN = nnz(NN(:,2));                                                 % Number of nodes with Heaviside enrichment
        IEN = nnz(NN(:,30));                                                % Number of inclusion nodes
        NEN = HEN+CTN+IEN;                                                  % Number of enriched nodes
        
        X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % X-coordinates of nodes
        Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Y-coordinates of nodes
        
        xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                    % Nodal coordinate matrix
        
        % Define the gauss points and the gauss weights
        if NEN == 4                                                         % Fully enriched element
            PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                     % Nodal crack level set values
            ZN = [ZETA(N1) ZETA(N2) ZETA(N3) ZETA(N4)];                     % Nodal inclusion level set values
            
            if HEN == 4
                [gp,gw,J] = subDomain(3,CN,PN,ZN,xyz,0,1,CCS);
            elseif CTN == 4
                [gp,gw,J] = subDomain(7,CN,PN,ZN,xyz,1,1,CCS);
            else
                [gp,gw,J] = subDomain(7,CN,PN,ZN,xyz,0,1,CCS);
            end
            
            if min(CN) < 0
                xi      = gp(:,1);
                eta     = gp(:,2);
                N       = 1/4*[(1-xi).*(1-eta) (1+xi).*(1-eta) (1+xi).*(1+eta) (1-xi).*(1+eta)];
                chi     = N(:,1)*CN(1)+N(:,2)*CN(2)+N(:,3)*CN(3)+N(:,4)*CN(4);
                R       = find(chi <= 0);
                J(R,:)  = [];
                gp(R,:) = [];
                gw(R,:) = [];
            end
            
        else                                                                % Partially enriched element
            if min(CN) < 0
                PN        = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];          % Nodal crack level set values
                ZN        = [ZETA(N1) ZETA(N2) ZETA(N3) ZETA(N4)];          % Nodal inclusion level set values
                
                [gp,gw,J] = subDomain(3,CN,PN,ZN,xyz,0,0,[]);
                xi        = gp(:,1);
                eta       = gp(:,2);
                N         = 1/4*[(1-xi).*(1-eta) (1+xi).*(1-eta) (1+xi).*(1+eta) (1-xi).*(1+eta)];
                chi       = N(:,1)*CN(1)+N(:,2)*CN(2)+N(:,3)*CN(3)+N(:,4)*CN(4);
                R         = find(chi <= 0);
                J(R,:)    = [];
                gp(R,:)   = [];
                gw(R,:)   = [];
            else
                [gp,gw] = gauss(6,'QUAD');
                J = [];
            end
        end
        
        % Loop through gauss points in current element to solve for J-integral
        Uenr = []; iGP = 1; iLoc = 1;
        for i = 1:length(gp)
            xi = gp(i,1); eta = gp(i,2);
            W  = gw(i);
            
            if isempty(J) == 0
                Ji   = [J(i,1) J(i,2);J(i,3) J(i,4)];                       % Jacobian for current subdivision
                detJ = det(Ji);                                             % Determinant of Jacobian for current subdivision
            else
                xyz1 = CCS*([X1 Y1]'-[xCT yCT]');                           % Crack tip coordinates of node 1
                xyz2 = CCS*([X2 Y2]'-[xCT yCT]');                           % Crack tip coordinates of node 2
                xyz3 = CCS*([X3 Y3]'-[xCT yCT]');                           % Crack tip coordinates of node 3
                xyz4 = CCS*([X4 Y4]'-[xCT yCT]');                           % Crack tip coordinates of node 4
                
                dxpdxi  = 1/4*(-(1-eta)*xyz1(1)+(1-eta)*xyz2(1)+(1+eta)*xyz3(1)-(1+eta)*xyz4(1));
                dxpdeta = 1/4*(-(1- xi)*xyz1(1)-(1+ xi)*xyz2(1)+(1+ xi)*xyz3(1)+(1- xi)*xyz4(1));
                dypdxi  = 1/4*(-(1-eta)*xyz1(2)+(1-eta)*xyz2(2)+(1+eta)*xyz3(2)-(1+eta)*xyz4(2));
                dypdeta = 1/4*(-(1- xi)*xyz1(2)-(1+ xi)*xyz2(2)+(1+ xi)*xyz3(2)+(1- xi)*xyz4(2));
                
                Je   = [dxpdxi dypdxi;dxpdeta dypdeta];                     % Jacobian for current element
                detJ = det(Je);                                             % Determinant of the Jacobian
            end
            
            % Define quadrilateral shape functions and derivatives
            N  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
            Nx = 2/lXElem*1/4*[-(1-eta);1-eta;1+eta;-(1+eta)];
            Ny = 2/lXElem*1/4*[-(1-xi);-(1+xi);1+xi;1-xi];
            
            Xgp = N(1)*X1+N(2)*X2+N(3)*X3+N(4)*X4;                          % The global X for the current gauss point
            Ygp = N(1)*Y1+N(2)*Y2+N(3)*Y3+N(4)*Y4;                          % The global Y for the current gauss point
            Zgp = N(1)*ZETA(N1)+N(2)*ZETA(N2)+N(3)*ZETA(N3)+N(4)*ZETA(N4);  % The value of zeta for the current gauss point
            
            X     = Xgp-xCT;                                                % Horizontal distance from crack tip to gauss point
            Y     = Ygp-yCT;                                                % Vertical distance from crack tip to gauss point
            XYloc = CCS*[X Y]';                                             % Change to crack tip coordinate system
            r     = sqrt(XYloc(1)^2+XYloc(2)^2);                            % Radius from crack tip to current gauss point
            theta = atan2(XYloc(2),XYloc(1));                               % Angle from crack tip to current gauss point
            
            U = [DISPLACEMENT(2*N1-1) DISPLACEMENT(2*N1) DISPLACEMENT(2*N2-1) DISPLACEMENT(2*N2)...
                DISPLACEMENT(2*N3-1) DISPLACEMENT(2*N3) DISPLACEMENT(2*N4-1) DISPLACEMENT(2*N4)];
            
            
            Benr = [];
            Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0   Nx(4)   0;...
                    0   Ny(1)   0   Ny(2)   0   Ny(3)   0   Ny(4);...
                  Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];
            
            iB = 1;
            for iN = 1:4
                if NN(iN,2) ~= 0                                            % Heaviside node
                    psi1 = PSI(N1);                                         % Value of phi at Node 1
                    psi2 = PSI(N2);                                         % Value of phi at Node 2
                    psi3 = PSI(N3);                                         % Value of phi at Node 3
                    psi4 = PSI(N4);                                         % Value of phi at Node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Calculate phi at the current gauss point
                    Hgp  = sign(psi);                                       % Find the value of H at the current gauss point
                    
                    Hi = NN(iN,3);                                          % Find the value of H at the current node
                    H  = Hgp-Hi;                                            % Find the value of H
                    
                    Ba = [Nx(iN)*H     0;
                              0    Ny(iN)*H;
                          Ny(iN)*H Nx(iN)*H];
                    Benr(:,iB:(iB+1)) = Ba;
                    iB = iB+2;
                    
                    if iGP == 1
                        Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NN(iN,2)-1) DISPLACEMENT(2*NN(iN,2))];
                        iLoc = iLoc+2;
                    end
                elseif NN(iN,12) ~= 0                                       % Bimaterial crack tip node
                    b = (Gm*(kf-1)-Gf*(km-1))/(Gm*(kf+1)+Gf*(km+1));        % Second Dundur's parameter
                    e = 1/(2*pi)*log((1-b)/(1+b));                          % Material constant
                    
                    % Common variables
                    sr  = sqrt(r);
                    st  = sin(theta);
                    st2 = sin(theta/2);
                    ct2 = cos(theta/2);
                    c   = 1/2/sqrt(r);
                    co  = CCS(1,1);
                    so  = CCS(1,2);
                    
                    a1gp  = sr*cos(e*log(r))*exp(-e*theta)*st2;             % Alpha 1 crack tip enrichment value
                    a2gp  = sr*cos(e*log(r))*exp(-e*theta)*ct2;             % Alpha 2 crack tip enrichment value
                    a3gp  = sr*cos(e*log(r))*exp(e*theta)*st2;              % Alpha 3 crack tip enrichment value
                    a4gp  = sr*cos(e*log(r))*exp(e*theta)*ct2;              % Alpha 4 crack tip enrichment value
                    a5gp  = sr*cos(e*log(r))*exp(e*theta)*st2*st;           % Alpha 5 crack tip enrichment value
                    a6gp  = sr*cos(e*log(r))*exp(e*theta)*ct2*st;           % Alpha 6 crack tip enrichment value
                    a7gp  = sr*sin(e*log(r))*exp(-e*theta)*st2;             % Alpha 7 crack tip enrichment value
                    a8gp  = sr*sin(e*log(r))*exp(-e*theta)*ct2;             % Alpha 8 crack tip enrichment value
                    a9gp  = sr*sin(e*log(r))*exp(e*theta)*st2;              % Alpha 9 crack tip enrichment value
                    a10gp = sr*sin(e*log(r))*exp(e*theta)*ct2;              % Alpha 10 crack tip enrichment value
                    a11gp = sr*sin(e*log(r))*exp(e*theta)*st2*st;           % Alpha 11 crack tip enrichment value
                    a12gp = sr*sin(e*log(r))*exp(e*theta)*ct2*st;           % Alpha 12 crack tip enrichment value
                    
                    a = [a1gp-NN(iN,5);...                                  % Shifted alpha 1 enrichment value
                         a2gp-NN(iN,7);...                                  % Shifted alpha 2 enrichment value
                         a3gp-NN(iN,9);...                                  % Shifted alpha 3 enrichment value
                         a4gp-NN(iN,11);...                                 % Shifted alpha 4 enrichment value
                         a5gp-NN(iN,13);...                                 % Shifted alpha 5 enrichment value
                         a6gp-NN(iN,15);...                                 % Shifted alpha 6 enrichment value
                         a7gp-NN(iN,17);...                                 % Shifted alpha 7 enrichment value
                         a8gp-NN(iN,19);...                                 % Shifted alpha 8 enrichment value
                         a9gp-NN(iN,21);...                                 % Shifted alpha 9 enrichment value
                         a10gp-NN(iN,23);...                                % Shifted alpha 10 enrichment value
                         a11gp-NN(iN,25);...                                % Shifted alpha 11 enrichment value
                         a12gp-NN(iN,27)];                                  % Shifted alpha 12 enrichment value
                    
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
                    Px = [px(1)*co+py(1)*-so;...
                          px(2)*co+py(2)*-so;...
                          px(3)*co+py(3)*-so;...
                          px(4)*co+py(4)*-so;...
                          px(5)*co+py(5)*-so;...
                          px(6)*co+py(6)*-so;...
                          px(7)*co+py(7)*-so;...
                          px(8)*co+py(8)*-so;...
                          px(9)*co+py(9)*-so;...
                          px(10)*co+py(10)*-so;...
                          px(11)*co+py(11)*-so;...
                          px(12)*co+py(12)*-so];
                    
                    % Derivative of bimaterial crack tip enrichment functions with respect to Y
                    Py = [px(1)*so+py(1)*co;...
                          px(2)*so+py(2)*co;...
                          px(3)*so+py(3)*co;...
                          px(4)*so+py(4)*co;...
                          px(5)*so+py(5)*co;...
                          px(6)*so+py(6)*co;...
                          px(7)*so+py(7)*co;...
                          px(8)*so+py(8)*co;...
                          px(9)*so+py(9)*co;...
                          px(10)*so+py(10)*co;...
                          px(11)*so+py(11)*co;...
                          px(12)*so+py(12)*co];
                    
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
                elseif NN(iN,4) ~= 0                                        % Crack tip node
                    a1gp = sqrt(r)*sin(theta/2);                            % Node 1 crack tip enrichment value
                    a2gp = sqrt(r)*cos(theta/2);                            % Node 2 crack tip enrichment value
                    a3gp = sqrt(r)*sin(theta)*sin(theta/2);                 % Node 3 crack tip enrichment value
                    a4gp = sqrt(r)*sin(theta)*cos(theta/2);                 % Node 4 crack tip enrichment value
                    
                    a1N = NN(iN,5); a2N = NN(iN,7); a3N = NN(iN,9); a4N = NN(iN,11);
                    a1  = a1gp-a1N; a2  = a2gp-a2N; a3  = a3gp-a3N; a4  = a4gp-a4N;
                    
                    c = 1/2/sqrt(r); ct = CCS(1,1); st = CCS(1,2);
                    
                    % Derivate of Phi_alpha with respect to X
                    Px = c*[-sin(theta/2)*ct              + cos(theta/2)*-st;...
                             cos(theta/2)*ct              + sin(theta/2)*-st;...
                            -sin(3*theta/2)*sin(theta)*ct + (sin(theta/2)+sin(3*theta/2)*cos(theta))*-st;...
                            -cos(3*theta/2)*sin(theta)*ct + (cos(theta/2)+cos(3*theta/2)*cos(theta))*-st];
                    
                    % Derivative of Phi_alpha with respect to Y
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
                end
            end
            
            if Zgp > 0, C = Cm; G = Gm; k = km; else C = Cf; G = Gf; k = kf; end
            
            B  = [Bu Benr];
            Xe = [U Uenr]';
            lB = size(B,2);
            
            % Solve for the stress and strain at the current gauss point
            strain = B*Xe;                                                  % Strain at current gauss point
            stress = C*strain;                                              % Stress at current gauss point
            
            % Derivates of q with respect to X and Y
            q     = qNode(iElem,:);                                         % q values at nodes of current element
            gradq = q*[Nx(1) Ny(1);Nx(2) Ny(2);Nx(3) Ny(3);Nx(4) Ny(4)];    % The derivative of q with respect to X and Y
            
            % Derivatives of nodal displacements with respect to X and Y
            Ux = B(1,1:2:lB)*Xe(1:2:lB);                                    % Derivative of U with respect to X
            Uy = B(2,2:2:lB)*Xe(1:2:lB);                                    % Derivative of U with respect to Y
            Vx = B(1,1:2:lB)*Xe(2:2:lB);                                    % Derivative of V with respect to X
            Vy = B(2,2:2:lB)*Xe(2:2:lB);                                    % Derivative of V with respect to Y
            
            % Convert quanties into crack tip coordinate system
            GradQ     = CCS*gradq';                                         % Gradient of q in crack tip coordinate system
            GradDisp  = CCS*[Ux Uy;Vx Vy]*CCS';                             % Gradient of displacement in crack tip coordinate system
            CalStress = CCS*[stress(1) stress(3);stress(3) stress(2)]*CCS'; % Stresses in crack tip coordinate system
            
            K1 = 1.0;                                                       % Chosen KI for pure Mode I asymptotic field
            K2 = 1.0;                                                       % Chosen KII for pure Mode II asymptotic field
            if nnz(NODES(:,12)) == 0                                        % Traditional crack auxiliary fields
                % Predefine variables to be used repeatidly
                SQR  = sqrt(r);                                             % The square-root of r
                CT   = cos(theta);                                          % The cosine of theta
                ST   = sin(theta);                                          % The sine of theta
                CT2  = cos(theta/2);                                        % The cosine of one half theta
                ST2  = sin(theta/2);                                        % The sine of one half theta
                C3T2 = cos(3*theta/2);                                      % The cosine of three halves theta
                S3T2 = sin(3*theta/2);                                      % The sine of three halves theta
                
                drdx =  CT;                                                 % Derivative of r with respect to X
                drdy =  ST;                                                 % Derivative of r with respect to Y
                dtdx = -ST/r;                                               % Derivative of theta with respect to X
                dtdy =  CT/r;                                               % Derivative of theta with respect to Y
                
                cAuxStress = sqrt(1/(2*pi));                                % Constant for auxiliary stress calculation
                cAuxDisp   = sqrt(1/(2*pi))/(2*G);                          % Constant for auxiliary displacement calculation
                
                AuxStress   = zeros(2,2);                                   % Initialize auxiliary stress matrix
                AuxGradDisp = zeros(2,2);                                   % Initialize gradient of displacement matrix
                for mode = 1:2
                    if mode == 1                                            % K1 = 1.0 and K2 = 0.0
                        AuxStress(1,1) = K1*cAuxStress/SQR*CT2*(1-ST2*S3T2);% Auxiliary stress component for Mode I loading in x-direction
                        AuxStress(2,2) = K1*cAuxStress/SQR*CT2*(1+ST2*S3T2);% Auxiliary stress component for Mode I loading in y-direction
                        AuxStress(1,2) = K1*cAuxStress/SQR*ST2*CT2*C3T2;    % Auxiliary shear stress component for Mode I loading
                        AuxStress(2,1) = AuxStress(1,2);                    % Auxiliary shear stress component for Mode I loading
                        
                        du1dr = K1*cAuxDisp*0.5/SQR*CT2*(k-CT);             % Derivative of auxiliary x-displacement with respect to r
                        du1dt = K1*cAuxDisp*SQR*(-0.5*ST2*(k-CT)+CT2*ST);   % Derivative of auxiliary x-displacement with respect to theta
                        du2dr = K1*cAuxDisp*0.5/SQR*ST2*(k-CT);             % Derivative of auxiliary y-displacement with respect to r
                        du2dt = K1*cAuxDisp*SQR*(0.5*CT2*(k-CT)+ST2*ST);    % Derivative of auxiliary y-displacement with respect to theta
                        
                        AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx;           % Auxiliary displacement gradient of u1 with respect to x
                        AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy;           % Auxiliary displacement gradient of u1 with respect to y
                        AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx;           % Auxiliary displacement gradient of u2 with respect to x
                        AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy;           % Auxiliary displacement gradient of u2 with respect to y
                        
                        AuxStrain(1,1) = AuxGradDisp(1,1);
                        AuxStrain(1,2) = 1/2*(AuxGradDisp(1,2)+AuxGradDisp(2,1));
                        AuxStrain(2,1) = AuxStrain(1,2);
                        AuxStrain(2,2) = AuxGradDisp(2,2);
                    elseif mode == 2                                         % K1 = 0.0 and K2 = 1.0
                        AuxStress(1,1) = -K2*cAuxStress/SQR*ST2*(2+CT2*C3T2);% Auxiliary stress component for Mode II loading in x-direction
                        AuxStress(2,2) =  K2*cAuxStress/SQR*ST2*CT2*C3T2;    % Auxiliary stress component for Mode II loading in y-direction
                        AuxStress(1,2) =  K2*cAuxStress/SQR*CT2*(1-ST2*S3T2);% Auxiliary shear stress component for Mode II loading
                        AuxStress(2,1) =  AuxStress(1,2);                    % Auxiliary shear stress component for Mode II loading
                        
                        du1dr =  K2*cAuxDisp*0.5/SQR*ST2*(k+2+CT);          % Derivative of auxiliary x-displacement with respect to r
                        du1dt =  K2*cAuxDisp*SQR*(0.5*CT2*(k+2+CT)-ST2*ST); % Derivative of auxiliary x-displacement with respect to theta
                        du2dr = -K2*cAuxDisp*0.5/SQR*CT2*(k-2+CT);          % Derivative of auxiliary y-displacement with respect to r
                        du2dt = -K2*cAuxDisp*SQR*(-0.5*ST2*(k-2+CT)-CT2*ST);% Derivative of auxiliary y-displacement with respect to theta
                        
                        AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx;           % Auxiliary displacement gradient of u1 with respect to x
                        AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy;           % Auxiliary displacement gradient of u1 with respect to y
                        AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx;           % Auxiliary displacement gradient of u2 with respect to x
                        AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy;           % Auxiliary displacement gradient of u2 with respect to y
                        
                        AuxStrain(1,1) = AuxGradDisp(1,1);
                        AuxStrain(1,2) = 1/2*(AuxGradDisp(1,2)+AuxGradDisp(2,1));
                        AuxStrain(2,1) = AuxStrain(1,2);
                        AuxStrain(2,2) = AuxGradDisp(2,2);
                    end
                    
                    I1 = (CalStress(1,1)*AuxGradDisp(1,1)+CalStress(2,1)*AuxGradDisp(2,1))*GradQ(1)+...
                        (CalStress(1,2)*AuxGradDisp(1,1)+CalStress(2,2)*AuxGradDisp(2,1))*GradQ(2);
                    
                    I2 = (AuxStress(1,1)*GradDisp(1,1)+AuxStress(2,1)*GradDisp(2,1))*GradQ(1)+...
                        (AuxStress(1,2)*GradDisp(1,1)+AuxStress(2,2)*GradDisp(2,1))*GradQ(2);
                    
                    if (isnan(I1) == 1) || (isnan(I2) == 1), break, end
                    
                    StrainEnergy = 0;
                    for j = 1:2
                        for k = 1:2
                            StrainEnergy = StrainEnergy+CalStress(j,k)*AuxStrain(j,k);
                        end
                    end
                    
                    I(mode,1) = I(mode,1)+(I1+I2-StrainEnergy*GradQ(1))*detJ*W;
                end
                iGP = iGP + 1;
            else                                                            % Bimaterial crack
                b = (Gf*(km-1)-Gm*(kf-1))/(Gf*(km+1)+Gm*(vf+1));            % Second Dundur's parameter
                e = 1/(2*pi)*log((1-b)/(1+b));                              % Material constant
                if Zgp > 0
                    Eo = Em; vo = vm; Go = Gm;
                elseif Zgp < 0
                    Eo = Ef; vo = vf; Go = Gf;
                end
                
                % Constants
                if ((theta >= 0) && (theta < pi)) || ((theta < -pi) && (theta >= -2*pi))
                    d = exp(-(pi-theta)*e);                                 % Constant, delta, upper half plane
                elseif ((theta >= pi) && (theta < 2*pi)) || ((theta < 0) && (theta >= -pi))
                    d = exp((pi+theta)*e);                                  % Constant, delta, lower half plane
                end
                
                BTA  = (0.5*cos(e*log(r))+e*sin(e*log(r)))/(0.25+e^2);      % Constant, beta
                BTAP = (0.5*sin(e*log(r))-e*cos(e*log(r)))/(0.25+e^2);      % Constant, beta prime
                GMA  = k*d-1/d;                                             % Constant, gamma
                GMAP = k*d+1/d;                                             % Constant, gamma prime
                PHI  = e*log(r)+theta/2;                                    % Constant, phi
                
                A = 1/(4*Go*cosh(e*pi));                                    % Constant, A
                B = sqrt(r/2/pi);                                           % Constant, B
                C = BTAP*GMA*cos(theta/2)-BTA*GMAP*sin(theta/2);            % Constant, C
                D = BTA*GMA*cos(theta/2)+BTAP*GMAP*sin(theta/2);            % Constant, D
                E = BTAP*GMAP*cos(theta/2)-BTA*GMA*sin(theta/2);            % Constant, E
                F = BTA*GMAP*cos(theta/2)+BTAP*GMA*sin(theta/2);            % Constant, F
                
                dCdr = e*D/r;                                               % Derivative of constant C with respect to r
                dCdt = -F/2+e*E;                                            % Derivative of constant C with respect to theta
                dDdr = -e*C/r;                                              % Derivative of constant D with respect to r
                dDdt = E/2+e*F;                                             % Derivative of constant D with respect to theta
                
                T1 = 2*d*sin(theta)*sin(PHI);                               % Constant, T1
                T2 = 2*d*sin(theta)*cos(PHI);                               % Constant, T2
                T3 = 2*d*cos(theta)*sin(PHI);                               % Constant, T3
                T4 = 2*d*cos(theta)*cos(PHI);                               % Constant, T4
                
                dT1dr = e*T2/r;                                             % Derivative of constant T1 with respect to r
                dT1dt = e*T1+T2/2+T3;                                       % Derivative of constant T1 with respect to theta
                dT2dr = -e*T1/r;                                            % Derivative of constant T2 with respect to r
                dT2dt = e*T2-T1/2+T4;                                       % Derivative of constant T2 with respect to theta
                
                drdx = cos(theta);                                          % Derivative of r with respect to x
                drdy = sin(theta);                                          % Derivative of r with respect to y
                dtdx = -sin(theta)/r;                                       % Derivative of theta with respect to x
                dtdy = cos(theta)/r;                                        % Derivative of theta with respect to y
                
                for mode = 1:2
                    if mode == 1                                            % K1 = 1.0 and K2 = 0.0
                        f1 = D+T1;                                          % Function, f1
                        f2 = -C-T2;                                         % Function, f2
                        
                        df1dr = dDdr+dT1dr;                                 % Derivative of function f1 with respect to r
                        df1dt = dDdt+dT1dt;                                 % Derivative of function f1 with respect to theta
                        df2dr = -dCdr-dT2dr;                                % Derivative of function f2 with respect to r
                        df2dt = -dCdt-dT2dt;                                % Derivative of function f2 with respect to theta
                        
                        df1dx = df1dr*drdx+df1dt*dtdx;                      % Derivative of function f1 with respect to x
                        df1dy = df1dr*drdy+df1dt*dtdy;                      % Derivative of function f1 with respect to y
                        df2dx = df2dr*drdx+df2dt*dtdx;                      % Derivative of function f2 with respect to x
                        df2dy = df2dt*drdy+df2dt*dtdy;                      % Derivative of function f2 with respect to y
                        
                        AuxGradDisp(1,1) = A*(B*df1dx+drdx*f1/(4*pi*B));    % Auxiliary displacement gradient of u1 with respect to x
                        AuxGradDisp(1,2) = A*(B*df1dy+drdy*f1/(4*pi*B));    % Auxiliary displacement gradient of u1 with respect to y
                        AuxGradDisp(2,1) = A*(B*df2dx+drdx*f2/(4*pi*B));    % Auxiliary displacement gradient of u2 with respect to x
                        AuxGradDisp(2,2) = A*(B*df2dy+drdy*f2/(4*pi*B));    % Auxiliary displacement gradient of u2 with respect to y
                        
                        AuxStrain(1,1) = AuxGradDisp(1,1);
                        AuxStrain(1,2) = 1/2*(AuxGradDisp(1,2)+AuxGradDisp(2,1));
                        AuxStrain(2,1) = AuxStrain(1,2);
                        AuxStrain(2,2) = AuxGradDisp(2,2);
                        
                        if plane == 1                                       % Plane stress
                            AuxStress(1,1) = Eo/(1-vo^2)*(AuxStrain(1,1)+vo*AuxStrain(2,2));
                            AuxStress(1,2) = Eo/(1+vo)*AuxStrain(1,2);
                            AuxStress(2,1) = Eo/(1+vo)*AuxStrain(2,1);
                            AuxStress(2,2) = Eo/(1-vo^2)*(vo*AuxStrain(1,1)+AuxStrain(2,2));
                        elseif plane == 2                                   % Plane strain
                            AuxStress(1,1) = Eo/(1+vo)/(1-2*vo)*((1-vo)*AuxStrain(1,1)+vo*AuxStrain(2,2));
                            AuxStress(1,2) = Eo/(1+vo)*AuxStrain(1,2);
                            AuxStress(2,1) = Eo/(1+vo)*AuxStrain(2,1);
                            AuxStress(2,2) = Eo/(1+vo)/(1-2*vo)*(vo*AuxStrain(1,1)+(1-vo)*AuxStrain(2,2));
                        end
                    elseif mode == 2                                        % K1 = 0.0 and K2 = 1.0
                        f1 = -C+T2;                                         % Function, f1
                        f2 = -D+T1;                                         % Function, f2
                        
                        df1dr = -dCdr+dT2dr;                                % Derivative of f1 with respect to r
                        df1dt = -dCdt+dT2dt;                                % Derivative of f1 with respect to theta
                        df2dr = -dDdr+dT1dr;                                % Derivative of f2 with respect to t
                        df2dt = -dDdt+dT1dt;                                % Derivative of f2 with respect to theta
                        
                        df1dx = df1dr*drdx+df1dt*dtdx;                      % Derivative of f1 with respect to x
                        df1dy = df1dr*drdy+df1dt*dtdy;                      % Derivative of f1 with respect to y
                        df2dx = df2dr*drdx+df2dt*dtdx;                      % Derivative of f2 with respect to x
                        df2dy = df2dt*drdy+df2dt*dtdy;                      % Derivative of f2 with respect to y
                        
                        AuxGradDisp(1,1) = A*(B*df1dx+drdx*f1/(4*pi*B));    % Auxiliary displacement gradient of u1 with respect to x
                        AuxGradDisp(1,2) = A*(B*df1dy+drdy*f1/(4*pi*B));    % Auxiliary displacement gradient of u1 with respect to y
                        AuxGradDisp(2,1) = A*(B*df2dx+drdx*f2/(4*pi*B));    % Auxiliary displacement gradient of u2 with respect to x
                        AuxGradDisp(2,2) = A*(B*df2dy+drdy*f2/(4*pi*B));    % Auxiliary displacement gradient of u2 with respect to y
                        
                        AuxStrain(1,1) = AuxGradDisp(1,1);
                        AuxStrain(1,2) = 1/2*(AuxGradDisp(1,2)+AuxGradDisp(2,1));
                        AuxStrain(2,1) = AuxStrain(1,2);
                        AuxStrain(2,2) = AuxGradDisp(2,2);
                        
                        if plane == 1                                       % Plane stress
                            AuxStress(1,1) = Eo/(1-vo^2)*(AuxStrain(1,1)+vo*AuxStrain(2,2));
                            AuxStress(1,2) = Eo/(1+vo)*AuxStrain(1,2);
                            AuxStress(2,1) = Eo/(1+vo)*AuxStrain(2,1);
                            AuxStress(2,2) = Eo/(1-vo^2)*(vo*AuxStrain(1,1)+AuxStrain(2,2));
                        elseif plane == 2                                   % Plane strain
                            AuxStress(1,1) = Eo/(1+vo)/(1-2*vo)*((1-vo)*AuxStrain(1,1)+vo*AuxStrain(2,2));
                            AuxStress(1,2) = Eo/(1+vo)*AuxStrain(1,2);
                            AuxStress(2,1) = Eo/(1+vo)*AuxStrain(2,1);
                            AuxStress(2,2) = Eo/(1+vo)/(1-2*vo)*(vo*AuxStrain(1,1)+(1-vo)*AuxStrain(2,2));
                        end
                    end
                    
                    I1 = (CalStress(1,1)*AuxGradDisp(1,1)+CalStress(2,1)*AuxGradDisp(2,1))*GradQ(1)+...
                         (CalStress(1,2)*AuxGradDisp(1,1)+CalStress(2,2)*AuxGradDisp(2,1))*GradQ(2);
                    
                    I2 = (AuxStress(1,1)*GradDisp(1,1)+AuxStress(2,1)*GradDisp(2,1))*GradQ(1)+...
                         (AuxStress(1,2)*GradDisp(1,1)+AuxStress(2,2)*GradDisp(2,1))*GradQ(2);
                    
                    StrainEnergy = 0;
                    for j = 1:2
                        for k = 1:2
                            StrainEnergy = StrainEnergy+CalStress(j,k)*AuxStrain(j,k);
                        end
                    end
                    
                    I(mode,1) = I(mode,1)+(I1+I2-StrainEnergy*GradQ(1))*detJ*W;
                end
                iGP = iGP + 1;
            end
        end
    end
    
    if isempty(setdiff(NODES(:,12),0)) == 1                                 % Traditional Crack
        % Find the effective modulus
        if plane == 1                                                       % Plane stress
            Eeff = Em;
        elseif plane == 2
            Eeff  = Em/(1-vm^2);                                            % Plane strain
        end
        
        % Solve for the mixed-mode stress intensity factors
        Kcalc   = I*Eeff/2;                                                 % Solve for Mode I and Mode II SIF
        KI(iJ)  = Kcalc(1);                                                 % Mode I  SIF
        KII(iJ) = Kcalc(2);                                                 % Mode II SIF
    else                                                                    % Bimaterial crack
        % Find the effective modulus
        if plane == 1                                                       % Plane stress
            Emeff = Em;
            Efeff = Ef;
        elseif plane == 2                                                   % Plane strain
            Emeff = Em/(1-vm^2);
            Efeff = Ef/(1-vf^2);
        end
        Eeff = 2*Emeff*Efeff/(Emeff+Efeff);
        
        % Solve for the mixed-mode stress intensity factors
        Kcalc   = I*Eeff/2*cosh(e*pi)*cosh(e*pi);                           % Solve for Mode I and Mode II SIF
        KI(iJ)  = Kcalc(1);                                                 % Mode I  SIF
        KII(iJ) = Kcalc(2);                                                 % Mode II SIF
    end
end