% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function plotStress(stressXX,stressXY,stressYY,stressVM)
% This function plots the stress for quadrilateral elements calculated in
% elemStress.m.

global CONNEC CRACK DOMAIN INC VOID XYZ

nElem  = DOMAIN(1)*DOMAIN(2);                                               % Number of elements
nXElem = DOMAIN(1);
nYElem = DOMAIN(2);
lXElem = DOMAIN(3);

for iPlot = 1:4
    if iPlot == 1
        Es = stressXX;                                                      % Plot sigmaXX
    elseif iPlot == 2
        Es = stressXY;                                                      % Plot sigmaXY
    elseif iPlot == 3
        Es = stressYY;                                                      % Plot sigmaYY
    elseif iPlot == 4
        Es = stressVM;                                                      % Plot sigmaVM
    end
    
    figure; hold on;
    
    smin = min(min(Es));                                                    % Find the minimum stress value
    smax = max(max(Es));                                                    % Find the maximum stress value

    if (smin < 1E-3) && (smax < 1E-3)                                       % Check if stresses are simply numerical noise
        if size(Es,2) == 1                                                  % Averaged stress values
            Es = zeros(length(Es));                                         % Define stress to be zero
        elseif size(Es,2) == 4                                              % Actual stress values
            Es = zeros(length(Es),4);                                       % Define stress to be zero
        end
        caxis([-0.1 0.1]);                                                  % Define bounds if stress is zero
    elseif (smin - smax) < 1E-6
        caxis([smin-0.1*smin smax+0.1*smax]);                               % Define bounds if stress is nonzero, constant
    else
        caxis([smin smax]);                                                 % Define bounds if stress is not constant
    end

    axis equal; axis off; colorbar('horiz');
    
    Ex = zeros(4,nElem); Ey = Ex; Ec = Ex;                                  % Initialize plotting matrices

    for iElem = 1:nElem
        NN          = CONNEC(iElem,2:5);                                    % Nodes of current element
        Ex(:,iElem) = XYZ(NN',2);                                           % X-coordinates of nodes
        Ey(:,iElem) = XYZ(NN',3);                                           % Y-coordinates of nodes
        if size(Es,2) == 1                                                  % Averaged stress values
            Ec(:,iElem) = Es(NN');                                          % Stress values at nodes
        elseif size(Es,2) == 4                                              % Actual stress values
            Ec(:,iElem) = Es(iElem,:)';                                     % Stress values at nodes
        end
    end
    
    % patch(Ex,Ey,Ec);
    H = patch(Ex,Ey,Ec);                                                    % Plot the stresses
    set(H,'LineStyle','none')
     
    if iPlot == 1
        title('\sigma_X_X');
    elseif iPlot == 2
        title('\sigma_X_Y');
    elseif iPlot == 3
        title('\sigma_Y_Y');
    elseif iPlot == 4
        title('\sigma_V_M');
    end

    % Plot the crack
    if isempty(CRACK) == 0
        nPt = size(CRACK,1);
        for iPt = 2:nPt
            x = [CRACK(iPt-1,1) CRACK(iPt,1)];
            y = [CRACK(iPt-1,2) CRACK(iPt,2)];
            plot(x,y,'w','LineWidth',2)
        end
    end
    
    % Plot the inclusions
    if isempty(INC) == 0
        nINC = size(INC,1);
        for i = 1:nINC
            xc = INC(i,1); yc = INC(i,2); rc = INC(i,3);
            theta = (0:256)*pi/2/256;
            plot(rc*cos(theta)+xc,rc*sin(theta)+yc,'w','LineWidth',2)
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
   
    hold off
end