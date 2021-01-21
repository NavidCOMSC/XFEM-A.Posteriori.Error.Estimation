function posteriori = FirstTryGetaposteriori(AvE, AvS, DUs, SS)

global DOMAIN CONNEC

nXElem  = DOMAIN(1);                                                        % Number of elements in the x-direction
nYElem  = DOMAIN(2);                                                        % Number of elements in the y-direction
lXElem  = DOMAIN(3);                                                        % Length of elements in the x-direction
DUs = DUs(:,[1,3,3,2]);
SS = SS(:,[1,3,3,2]);
AvE = AvE(:,[1,3,3,2]);
AvS = AvS(:,[1,3,3,2]);

%size_coordinates = (nXElem+1)*(nYElem+1);
posteriori = zeros(nXElem*nYElem,1);
area = (lXElem)^2;

for iElem = 1:(nXElem*nYElem)
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN = [N1 N2 N3 N4];
    
    for i=1:4
        %size(SS(iElem,i))
        %size(DUs(iElem,i))
        %size(AvS(NN,i)')
        %size(AvE(NN,i))
%         j=1;
%         eta4(j) = eta4(j) + Area4 * (Sigma4(j,i)*Eps4(j,i) ...
% 	+ AvS(elements4(j,:),i)'*[4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4]* ...
% 	AvE(elements4(j,:),i)/36 ...
% 	- AvS(elements4(j,:),i)'*[1;1;1;1]*Eps4(j,i)/4 ...
% 	- AvE(elements4(j,:),i)'*[1;1;1;1]*Sigma4(j,i)/4);
        %SS(iElem,:)
        %DUs(iElem,:)
        %AvS(NN,:)
        %AvE(NN,:)
        
        posteriori(iElem) = posteriori(iElem)+area*(SS(iElem,i)*DUs(iElem,i)...
            +AvS(NN,i)'*[4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4]*AvE(NN,i)/36 ... 
             - AvS(NN,i)'*[1;1;1;1]*DUs(iElem,i)/4 -...
                 AvE(NN,i)'*[1;1;1;1]*SS(iElem,i)/4);
             
             
    end
    %posteriori(iElem)
    %pause
end




