% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function plotContour(Sxx,Sxy,Syy,Svm)
% This file plots the stress contours from the stress values calculated by
% elemStress.m

global CONNEC CRACK DOMAIN INC PLOT VOID XYZ

nXElem = DOMAIN(1);                                                         % Number of elements in x-direction
nYElem = DOMAIN(2);                                                         % Number of elements in y-direction
nNodes = (nXElem+1)*(nYElem+1);                                             % Number of nodes
lXElem = DOMAIN(3);                                                         % Length of elements in x-direction

X = zeros(nXElem+1,nYElem+1); Y = X; Z = X;

% Average nodal values if not already done
if size(Sxx,2) == 4                                                         % Average nodal stress values
    stressXX = Sxx; stressXY = Sxy; stressYY = Syy; stressVM = Svm;
    Sxx = zeros(nNodes,2); Sxy = Sxx; Syy = Sxx; Svm = Sxx;
    
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
end

% Plot the contours
if PLOT(5,3) == 1, nPlot = 1; else nPlot = 3; end
    
for i = 1:nPlot
    figure; hold on;
    if nPlot == 3
        if i == 1
            z = Sxx;
            title('\sigma_x_x')
        elseif i == 2
            z = Sxy;
            title('\sigma_x_y')
        elseif i == 3
            z = Syy;
            title('\sigma_y_y')
        end
    else
        z = Svm;
        title('\sigma_v_m')
    end
    
    iNode = 1;
    for yE = 1:(nYElem+1)
        for xE = 1:(nXElem+1)
            X(xE,yE) = XYZ(iNode,2);
            Y(xE,yE) = XYZ(iNode,3);
            Z(xE,yE) = z(iNode);
            iNode = iNode+1;
        end
    end
    
    if PLOT(5,2) == 1                                                       % Filled contour
        contourf(X,Y,Z,40,'LineStyle','none');
    else
        contour(X,Y,Z);                                                     % Contour lines
    end
    
    % Plot the crack
    if isempty(CRACK) == 0
        nPt = size(CRACK,1);
        for iPt = 2:nPt
            xc = [CRACK(iPt-1,1) CRACK(iPt,1)];
            yc = [CRACK(iPt-1,2) CRACK(iPt,2)];
            plot(xc,yc,'w','LineWidth',2)
        end
    end
    
    % Plot the inclusions
    if isempty(INC) == 0
        nINC = size(INC,1);
        xMax = nXElem*lXElem;
        yMax = nYElem*lXElem;
        for iI = 1:nINC
            if size(INC,2)==3
                xc = INC(iI,1); yc = INC(iI,2); rc = INC(iI,3);
                theta = (0:256)*pi*2/256;
                xp = rc*cos(theta)+xc;
                yp = rc*sin(theta)+yc;
                for j = 1:length(xp)
                    if (xp(j) > xMax) || (xp(j) < 0)
                        xp(j) = NaN; yp(j) = NaN;
                    elseif (yp(j) > yMax) || (yp(j) < 0)
                        xp(j) = NaN; yp(j) = NaN;
                    end
                end
                plot(xp,yp,'w','LineWidth',2)
            else
                line([INC(iI,1), INC(iI,3)], [INC(iI,2), INC(iI,4)], 'LineWidth', 2)
            end
                
            
        end
    end
    
    % Plot the voids
    if isempty(VOID) == 0
        nVOID = size(VOID,1);
        xMax = nXElem*lXElem;
        yMax = nYElem*lXElem;
        for iI = 1:nVOID
            xc = VOID(iI,1); yc = VOID(iI,2); rc = VOID(iI,3);
            theta = 0:0.5:360;
            xp = rc*cosd(theta)+xc;
            yp = rc*sind(theta)+yc;
            for j = 1:length(xp)
                if (xp(j) > xMax) || (xp(j) < 0)
                    xp(j) = 0; yp(j) = 0;
                elseif (yp(j) > yMax) || (yp(j) < 0)
                    xp(j) = 0; yp(j) = 0;
                end
            end
            fill(xp,yp,'w','LineStyle','none')
        end
    end
    
    axis equal; axis([0 nXElem*lXElem 0 nYElem*lXElem]);
    set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
    colorbar('horiz'); 
end
