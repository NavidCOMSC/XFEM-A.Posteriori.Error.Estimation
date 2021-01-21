% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function ctipNodes
% This function assigns nodes enriched with the crack tip enrichment
% function values associated with the near tip asymptotic displacement
% field.

global CRACK NODES PHI XYZ

iSeg = size(CRACK,1);                                                       % Number of crack segments
nCT  = size(PHI,2);                                                         % Number of crack tips

% Define coordinates of crack tip(s)
if nCT == 1
    xCT = CRACK(iSeg,1);                                                    % X-coordinate of crack tip
    yCT = CRACK(iSeg,2);                                                    % Y-coordinate of crack tip
elseif nCT == 2
    xCT = [CRACK(iSeg,1) CRACK(1,1)];                                       % X-coordinate of crack tips
    yCT = [CRACK(iSeg,2) CRACK(1,2)];                                       % Y-coordinate of crack tips
end

for iNode = 1:size(NODES,1)
    if NODES(iNode,4) ~= 0
        XN = XYZ(iNode,2);                                                  % Nodal x-coordinate
        YN = XYZ(iNode,3);                                                  % Nodal y-coordinate
        X  = XN-xCT;                                                        % Horizontal distance from crack tip(s) to current node
        Y  = YN-yCT;                                                        % Vertical distance from crack tip(s) to current node
        for i = 1:length(X)
            rt = sqrt(X(i)^2+Y(i)^2);
            tt = atan2(Y(i),X(i));
            if i == 1; r = rt; theta = tt; end
            if i > 1; if rt < r, r = rt; theta = tt; end, end
        end
        
        NODES(iNode,5)  = sqrt(r)*sin(theta/2);                             % Alpha 1 crack tip enrichment value
        NODES(iNode,7)  = sqrt(r)*cos(theta/2);                             % Alpha 2 crack tip enrichment value
        NODES(iNode,9)  = sqrt(r)*sin(theta)*sin(theta/2);                  % Alpha 3 crack tip enrichment value
        NODES(iNode,11) = sqrt(r)*sin(theta)*cos(theta/2);                  % Alpha 4 crack tip enrichment value
    end
end
