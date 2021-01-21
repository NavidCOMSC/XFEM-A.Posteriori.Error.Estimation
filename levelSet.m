% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [omega] = levelSet(iter)
% This function creates the ZETA level set function representing the
% materials within the domain. This function creates the CHI level set
% function representing voids within the domain. This function creates the
% PHI and PSI level set functions needed to represent the discontinuity
% created by the crack. The elements with Heaviside and crack tip
% enrichments are determined and node numbers are assigned to the nodes
% which have enriched degrees of freedom.

global CHI CONNEC CRACK DOMAIN INC NODES PHI PLOT PSI XYZ VOID ZETA

nXElem = DOMAIN(1);                                                         % Number of elements in the x-direction
nYElem = DOMAIN(2);                                                         % Number of elements in the y-direction
lXElem = DOMAIN(3);                                                         % Length of elements in the x-direction
nElem  = nXElem*nYElem;                                                     % Number of elements in domain
nNode  = (nXElem+1)*(nYElem+1);                                             % Number of nodes in domain
nPt    = size(CRACK,1);                                                     % Number of data points defining crack

if nnz(PLOT(1,7)) == 0
    nBand = (sqrt(2)+0.005)*lXElem;                                         % Radius of narrow band about crack
else
    nBand = PLOT(1,7)*lXElem;                                               % User-defined narrow band radius
end

if iter == 1, incNode = 1; else incNode = 1+max(NODES(:,2))-nNode; end

if iter == 1
    % Create the material level set function (ZETA)
    if isempty(INC) == 0
        if size(INC,2) == 3                                                 % Circular inclusions
            disp('Circular inclusion detected')
            nInc = size(INC,1);                                             % Number of inclusions
            ZETA = zeros(nInc,nNode);
            for iInc = 1:nInc
                xc = INC(iInc,1);                                           % x-coordinate of fiber center
                yc = INC(iInc,2);                                           % y-coordinate of fiber center
                rc = INC(iInc,3);                                           % radius of fiber
                
                % Define the level set for the inclusion
                for iNode = 1:nNode
                    xi = XYZ(iNode,2);                                      % X-coordinate of current node
                    yi = XYZ(iNode,3);                                      % Y-coordinate of current node
                    Z  = sqrt((xi-xc)^2+(yi-yc)^2)-rc;                      % Level set value corresponding to current inclusion
                    if abs(Z) < 1e-6
                        ZETA(iInc,iNode) = 0;
                    else
                        ZETA(iInc,iNode) = Z;
                    end
                end
            end
        else                                                                % Linear inclusion
            x0 = INC(1,1);                                                  % x-coordinate of first defining point
            y0 = INC(1,2);                                                  % y-coordinate of first defining point
            x1 = INC(1,3);                                                  % x-coordinate of second defining point
            y1 = INC(1,4);                                                  % y-coordinate of second defining point
            l  = sqrt((x1-x0)^2+(y1-y0)^2);                                 % Length of defining segment
            for iNode = 1:nNode                                             % Build zeta level set function at each node
                x = XYZ(iNode,2);                                           % X-coordinate of current node
                y = XYZ(iNode,3);                                           % Y-coordinate of current node
                c = (y0-y1)*x+(x1-x0)*y+(x0*y1-x1*y0);                      % Distance from node to zero level set of zeta
                ZETA(iNode) = c/l;                                          % Zeta level set value
            end
        end
        
        % Combine separate level set functions into a single global function
        ZETA = min(ZETA,[],1);
        
        % Find the elements which need the inclusion enrichment
        for iElem = 1:nElem
            Z1 = ZETA(CONNEC(iElem,2));
            Z2 = ZETA(CONNEC(iElem,3));
            Z3 = ZETA(CONNEC(iElem,4));
            Z4 = ZETA(CONNEC(iElem,5));
            Z  = [Z1 Z2 Z3 Z4];
            
            if max(Z)*min(Z) < 0
                for iNode = 2:5
                    if NODES(CONNEC(iElem,iNode),30) == 0
                        NODES(CONNEC(iElem,iNode),30) = nNode+incNode;
                        NODES(CONNEC(iElem,iNode),31) = Z(iNode-1);
                        incNode = incNode+1;
                    end
                end
            end
        end
    else
        ZETA  = ones(1,nNode);
    end
    
    % Create the level set function defining the voids
    if isempty(VOID) == 0
        nVoid = size(VOID,1);
        CHI   = zeros(nVoid,nNode);
        for iVoid = 1:nVoid
            xc = VOID(iVoid,1);                                             % x-coordinate of fiber center
            yc = VOID(iVoid,2);                                             % y-coordinate of fiber center
            rc = VOID(iVoid,3);                                             % radius of fiber
            
            % Define the level set for the inclusion
            for iNode = 1:nNode
                xi = XYZ(iNode,2);                                          % X-coordinate of current node
                yi = XYZ(iNode,3);                                          % Y-coordinate of current node
                C  = sqrt((xi-xc)^2+(yi-yc)^2)-rc;                          % Level set value corresponding to current inclusion
                if abs(C) < 1e-6
                    CHI(iVoid,iNode) = 0;
                else
                    CHI(iVoid,iNode) = C;
                end
            end
        end
        
        % Combine separate level set functions into a single global function
        CHI = min(CHI,[],1);
    else
        CHI = zeros(1,nNode);
    end
end

% Create the level set functions (PHI and PSI) defining the crack
if isempty(CRACK) == 0
    nCT = 2;                                                                % Default number of crack tips
    if     (CRACK(1,1)   == 0) || (CRACK(1,1)   == nXElem*lXElem)           % Check for edge crack
        nCT = nCT-1;
    elseif (CRACK(nPt,1) == 0) || (CRACK(nPt,1) == nXElem*lXElem)
        nCT = nCT-1;
    elseif (CRACK(1,2)   == 0) || (CRACK(1,2)   == nYElem*lXElem)
        nCT = nCT-1;
    elseif (CRACK(nPt,2) == 0) || (CRACK(nPt,2) == nYElem*lXElem)
        nCT = nCT-1;
    end
    
    if nCT == 2
        disc  = CRACK(nPt,:)-CRACK(nPt-1,:);                                % Horizontal and vertical distances for current crack segment
        omega = atan2(disc(2),disc(1));                                     % Crack angle with respect to horizontal
        disc  = CRACK(1,:)-CRACK(2,:);                                      % Horizontal and vertical distances for current crack segment
        omega(2) = atan2(disc(2),disc(1));                                  % Crack angle with respect to horizontal
    elseif nCT == 1
        disc  = CRACK(nPt,:)-CRACK(nPt-1,:);                                % Horizontal and vertical distances for current crack segment
        omega = atan2(disc(2),disc(1));                                     % Crack angle with respect to horizontal
    end
    
    % Find the elements to search for new level set values
    PHI = sparse(nNode,nCT);                                                % Initialize phi
    if iter == 1, PSI = sparse(nNode,1); end                                % Initialize psi
    
    for iCT = 1:nCT
        
        if iCT == 1
            dFSeg  = sqrt((CRACK(nPt-1,1)-CRACK(nPt,1))^2+(CRACK(nPt-1,2)-CRACK(nPt,2))^2);
            radius = dFSeg+nBand;                                               % Define radius of nodal search for defining level set
            xCT    = CRACK(nPt,1);                                              % X-coordinate of crack tip
            yCT    = CRACK(nPt,2);                                              % Y-coordinate of crack tip
            CCS    = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];% Convert from global to crack tip coordinate system
        elseif iCT == 2
            dFSeg  = sqrt((CRACK(1,1)-CRACK(2,1))^2+(CRACK(1,2)-CRACK(2,2))^2);
            radius = dFSeg+nBand;                                               % Define radius of nodal search for defining level set
            xCT    = CRACK(1,1);                                                % X-coordinate of crack tip
            yCT    = CRACK(1,2);                                                % Y-coordinate of crack tip
            CCS    = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];% Convert from global to crack tip coordinate system
        end
        
        dist = zeros(1,nNode);                                              % Initialize distance vector
        for iN = 1:nNode
            Xn       = XYZ(iN,2);                                           % X-coordinate for the current node
            Yn       = XYZ(iN,3);                                           % Y-coordinate for the current node
            X        = Xn-xCT;                                              % Horizontal distance from crack tip to current node
            Y        = Yn-yCT;                                              % Vertical distance from crack tip to current node
            XYloc    = CCS*[X Y]';                                          % Change to crack tip coordinates
            r        = sqrt(XYloc(1)^2+XYloc(2)^2);                         % Radius from crack tip to current gauss point
            dist(iN) = r;                                                   % Store radius value
        end
        
        temp   = dist-radius;                                               % Determine whether or not the node is outside the search radius
        domain = find(temp <= 0);                                           % Find nodes within search radius
        
        % Compute the PHI level set functions for the main crack tip(s)
        for iNode = 1:length(domain)
            cNode = domain(iNode);                                          % Current node within search radius
            x     = XYZ(cNode,2);                                           % X-coordinate for the current node
            y     = XYZ(cNode,3);                                           % Y-coordinate for the current node
            
            % Define phi for first crack tip
            disc = CRACK(nPt,:)-CRACK(nPt-1,:);                             % Horizontal and vertical distances for current crack segment
            t    = 1/norm(disc)*disc;                                       % Tangent to current crack segment
            phi  = ([x y]-CRACK(nPt,:))*t';
            if phi == 0, phi = 1e-6; end
            
            if nCT == 2
                disc = CRACK(1,:)-CRACK(2,:);                               % Horizontal and vertical distances for current crack segment
                t    = 1/norm(disc)*disc;                                   % Tangent to current crack segment
                phi(2)  = ([x y]-CRACK(1,:))*t';
                if phi(2) == 0, phi(2) = 1e-6; end
            end
            
            % Define phi and psi at nodes within narrow band
            if iCT == 1
                if abs(phi(1)) < nBand                                      % Check if phi is within narrow band
                    
                    % Define psi
                    dist = zeros(nPt-1,1); sine = dist;
                    for iSeg = 1:(nPt-1)
                        x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                        x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
                        
                        m = (y2-y1)/(x2-x1);
                        if isinf(m) == 1, m = 1e6; end
                        b = y1-m*x1;
                        
                        xo = (m*y+x-m*b)/(m^2+1);
                        yo = (m^2*y+m*x+b)/(m^2+1);
                        
                        if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                        if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                        
                        dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                        sine(iSeg) = sign(x2-x1)*sign(y-yo);
                    end
                    
                    dMin = min(abs(dist));
                    ind  = find(abs(dist) == dMin);
                    if length(ind) == 2, ind(2) = []; end
                    psi = dist(ind)*sine(ind);
                    if psi == 0, psi = 1e-6; end
                    
                    % Check if psi is within defined narrow band
                    if abs(psi) < nBand
                        PHI(cNode,1) = phi(1);
                        PSI(cNode,1) = psi;
                    end
                elseif phi(1) < 0                                           % Check if phi is negative
                    
                    % Define psi
                    dist = zeros(nPt-1,1); sine = dist;
                    for iSeg = 1:(nPt-1)
                        x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                        x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
                        
                        m = (y2-y1)/(x2-x1);
                        if isinf(m) == 1, m = 1e6; end
                        b = y1-m*x1;
                        
                        xo = (m*y+x-m*b)/(m^2+1);
                        yo = (m^2*y+m*x+b)/(m^2+1);
                        
                        if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                        if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                        
                        dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                        sine(iSeg) = sign(x2-x1)*sign(y-yo);
                    end
                    
                    dMin = min(abs(dist));
                    ind  = find(abs(dist) == dMin);
                    if length(ind) == 2, ind(2) = []; end
                    psi = dist(ind)*sine(ind);
                    if psi == 0, psi = 1e-6; end
                    
                    % Check if psi is within narrow band
                    if abs(psi) < nBand
                        PSI(cNode,1) = psi;
                    end
                end
            end
            
            % Define phi for second crack tip
            if iCT == 2
                if abs(phi(2)) < nBand                                      % Check if phi is within narrow band

                    % Define psi
                    dist = zeros(nPt-1,1); sine = dist;
                    for iSeg = 1:(nPt-1)
                        x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                        x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
                        
                        m = (y2-y1)/(x2-x1);
                        if isinf(m) == 1, m = 1e6; end
                        b = y1-m*x1;
                        
                        xo = (m*y+x-m*b)/(m^2+1);
                        yo = (m^2*y+m*x+b)/(m^2+1);
                        
                        if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                        if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                        
                        dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                        sine(iSeg) = sign(x2-x1)*sign(y-yo);
                    end
                    
                    dMin = min(abs(dist));
                    ind  = find(abs(dist) == dMin);
                    if length(ind) == 2, ind(2) = []; end
                    psi = dist(ind)*sine(ind);
                    if psi == 0, psi = 1e-6; end
                    
                    % Check if psi is within defined narrow band
                    if abs(psi) < nBand
                        if PHI(cNode,1) ~= 0
                            if abs(PHI(cNode,1)) > abs(phi(2))
                                PHI(cNode,1) = phi(2);
                                PSI(cNode,1) = psi;
                            end
                        else
                            PHI(cNode,1) = phi(2);
                            PSI(cNode,1) = psi;
                        end
                    end
                elseif phi(2) < 0                                           % Check if phi is negative
                    if phi(1) < 0
                        % Define psi
                        dist = zeros(nPt-1,1); sine = dist;
                        for iSeg = 1:(nPt-1)
                            x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                            x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
                            
                            m = (y2-y1)/(x2-x1);
                            if isinf(m) == 1, m = 1e6; end
                            b = y1-m*x1;
                            
                            xo = (m*y+x-m*b)/(m^2+1);
                            yo = (m^2*y+m*x+b)/(m^2+1);
                            
                            if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                            if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                            
                            dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                            sine(iSeg) = sign(x2-x1)*sign(y-yo);
                        end
                        
                        dMin = min(abs(dist));
                        ind  = find(abs(dist) == dMin);
                        if length(ind) == 2, ind(2) = []; end
                        psi = dist(ind)*sine(ind);
                        if psi == 0, psi = 1e-6; end
                        
                        % Check if psi is within narrow band
                        if abs(psi) < nBand
                            PSI(cNode,1) = psi;
                        end
                    end
                end
            end
        end
    end
    
    % Define the crack tip enriched nodes
    ctNodes  = [];
    bmNodes  = [];
    I        = find(PHI ~= 0);                                              % Nodes with defined phi
    [a,b,c1] = intersect(I,CONNEC(:,2)');                                   % Find elements with defined phi
    [a,b,c2] = intersect(I,CONNEC(:,3)');                                   % Find elements with defined phi
    [a,b,c3] = intersect(I,CONNEC(:,4)');                                   % Find elements with defined phi
    [a,b,c4] = intersect(I,CONNEC(:,5)');                                   % Find elements with defined phi
    ctElem   = unique([c1 c2 c3 c4]);                                       % Candidate elements for crack tip enrichment
    
    for iElem = 1:length(ctElem)
        cElem = ctElem(iElem);
        phiE  = PHI(CONNEC(cElem,2:5));
        if nnz(phiE) == 4
            psiE = PSI(CONNEC(cElem,2:5));
            if max(psiE)*min(psiE) < 1e-6
                if max(phiE)*min(phiE) < 1e-6
                    for iN = 1:4
                        ctNodes = [ctNodes CONNEC(cElem,iN+1)];
                    end
                end
            end
        end
    end
    NODES(ctNodes,4) = NaN;
    
    % Define Heaviside enriched nodes
    hNodes   = [];
    I        = find(PSI ~= 0);                                              % Nodes with defined psi
    [a,b,c1] = intersect(I,CONNEC(:,2)');                                   % Find elements with defined psi
    [a,b,c2] = intersect(I,CONNEC(:,3)');                                   % Find elements with defined psi
    [a,b,c3] = intersect(I,CONNEC(:,4)');                                   % Find elements with defined psi
    [a,b,c4] = intersect(I,CONNEC(:,5)');                                   % Find elements with defined psi
    hElem    = unique([c1 c2 c3 c4]);                                       % Candidate elements for Heaviside enrichment
    
    for iElem = 1:length(hElem)
        cElem = hElem(iElem);
        psiE  = PSI(CONNEC(cElem,2:5));
        if nnz(psiE) == 4
            if (max(psiE) == 1e-6) || (min(psiE) == 1e-6)
                for iN = 1:4
                    gN = CONNEC(cElem,iN+1);
                    if NODES(gN,4) == 0
                        if psiE(iN) == 1e-6
                            if PHI(gN) <= 0
                                if NODES(gN,2) == 0
                                    hNodes = [hNodes gN];
                                end
                            end
                        end
                    end
                end
            elseif max(psiE)*min(psiE) < 0
                for iN = 1:4
                    gN = CONNEC(cElem,iN+1);
                    if PHI(gN) <= 0
                        if NODES(gN,4) == 0
                            if NODES(gN,2) == 0
                                hNodes = [hNodes gN];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Preliminaries to numbering enriched nodes
    nNode   = nNode+1+incNode-1;
    ctNodes = unique(ctNodes);
    bmNodes = unique(bmNodes);
    hNodes  = unique(hNodes);
    
    % Number the Heaviside nodes
    for i = 1:length(hNodes)
        NODES(hNodes(i),2) = nNode;
        nNode = nNode+1;
    end
    
    % Number the crack tip nodes
    for i = 1:length(ctNodes)
        NODES(ctNodes(i),4)  = nNode;
        NODES(ctNodes(i),6)  = nNode+1;
        NODES(ctNodes(i),8)  = nNode+2;
        NODES(ctNodes(i),10) = nNode+3;
        nNode = nNode+4;
    end
    
    % Number the bimaterial crack tip nodes
    for i = 1:length(bmNodes)
        NODES(bmNodes(i),4)  = nNode;
        NODES(bmNodes(i),6)  = nNode+1;
        NODES(bmNodes(i),8)  = nNode+2;
        NODES(bmNodes(i),10) = nNode+3;
        NODES(bmNodes(i),12) = nNode+4;
        NODES(bmNodes(i),14) = nNode+5;
        NODES(bmNodes(i),16) = nNode+6;
        NODES(bmNodes(i),18) = nNode+7;
        NODES(bmNodes(i),20) = nNode+8;
        NODES(bmNodes(i),22) = nNode+9;
        NODES(bmNodes(i),24) = nNode+10;
        NODES(bmNodes(i),26) = nNode+11;
        nNode = nNode+12;
    end
else
    omega  = [];
end