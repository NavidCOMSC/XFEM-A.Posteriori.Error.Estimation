% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [Q,W] = gauss(npt,type)
% This function contains various gauss integration points and weights for
% quadrilateral and triangular elements.

if strcmp(type,'QUAD') == 1
    switch npt
        
        case 1
            np(1) = 0;
            nw(1) = 2;
        
        case 2
            np(1) = -0.577350269189626;
            np(2) =  0.577350269189626;

            nw(1) =  1.000000000000000;
            nw(2) =  1.000000000000000;
        case 3
            np(1) = -0.774596669241483;
            np(2) =  0.000000000000000;
            np(3) =  0.774596669241483;

            nw(1) =  0.555555555555556;
            nw(2) =  0.888888888888889;
            nw(3) =  0.555555555555556;
        case 4
            np(1) = -0.861136311594053;
            np(2) = -0.339981043584856;
            np(3) =  0.339981043584856;
            np(4) =  0.861136311594053;

            nw(1) =  0.347854845137454;
            nw(2) =  0.652145154862546;
            nw(3) =  0.652145154862546;
            nw(4) =  0.347854845137454;
        case 5
            np(1) = -0.906179845938664;
            np(2) = -0.538469310105683;
            np(3) =  0.000000000000000;
            np(4) =  0.538469310105683;
            np(5) =  0.906179845938664;

            nw(1) =  0.236926885056189;
            nw(2) =  0.478628670499366;
            nw(3) =  0.568888888888889;
            nw(4) =  0.478628670499366;
            nw(5) =  0.236926885056189;
        case 6
            np(1) = -0.932469514203152;
            np(2) = -0.661209386466265;
            np(3) = -0.238619186083197;
            np(4) =  0.238619186083197;
            np(5) =  0.661209386466265;
            np(6) =  0.932469514203152;

            nw(1) =  0.171324492379170;
            nw(2) =  0.360761573048139;
            nw(3) =  0.467913934572691;
            nw(4) =  0.467913934572691;
            nw(5) =  0.360761573048139;
            nw(6) =  0.171324492379170;
        case 7
            np(1) = -0.949107912342759;
            np(2) = -0.741531185599394;
            np(3) = -0.405845151377397;
            np(4) =  0.000000000000000;
            np(5) =  0.405845151377397;
            np(6) =  0.741531185599394;
            np(7) =  0.949109712342759;

            nw(1) =  0.129484966168870;
            nw(2) =  0.279705391489277;
            nw(3) =  0.381830050505119;
            nw(4) =  0.417959183673469;
            nw(5) =  0.381830050505119;
            nw(6) =  0.279705391489277;
            nw(7) =  0.129484966168870;
        case 8
            np(1) = -0.960289856497536;
            np(2) = -0.796666477413627;
            np(3) = -0.525532409916329;
            np(4) = -0.183434642495650;
            np(5) =  0.183434642495650;
            np(6) =  0.525532409916329;
            np(7) =  0.796666477413627;
            np(8) =  0.960289856497536;

            nw(1) =  0.101228536290376;
            nw(2) =  0.222381034453374;
            nw(3) =  0.313706645877887;
            nw(4) =  0.362683783378362;
            nw(5) =  0.362683783378362;
            nw(6) =  0.313706645877887;
            nw(7) =  0.222381034453374;
            nw(8) =  0.101228536290376;
        case 9
            np(1) = -0.968160239507626;
            np(2) = -0.836031107326636;
            np(3) = -0.613371432700590;
            np(4) = -0.324253423403809;
            np(5) =  0.000000000000000;
            np(6) =  0.324253423403809;
            np(7) =  0.613371432700590;
            np(8) =  0.836031107326636;
            np(9) =  0.968160239507626;

            nw(1) =  0.081274388361574;
            nw(2) =  0.180648160694857;
            nw(3) =  0.260610696402935;
            nw(4) =  0.312347077040003;
            nw(5) =  0.330239355001260;
            nw(6) =  0.312347077040003;
            nw(7) =  0.261610696402935;
            nw(8) =  0.180648160694857;
            nw(9) =  0.081274388361574;
        case 10
            np(1)  = -0.973906528517172;
            np(2)  = -0.865063366688985;
            np(3)  = -0.679409568299024;
            np(4)  = -0.433395394129247;
            np(5)  = -0.148874338981631;
            np(6)  =  0.148874338981631;
            np(7)  =  0.433395394129247;
            np(8)  =  0.679409568299024;
            np(9)  =  0.865063366688985;
            np(10) =  0.973906528517172;

            nw(1)  =  0.066671344308688;
            nw(2)  =  0.149451349150581;
            nw(3)  =  0.219086362515982;
            nw(4)  =  0.269266719309996;
            nw(5)  =  0.295524224714753;
            nw(6)  =  0.295524224714753;
            nw(7)  =  0.269266719309996;
            nw(8)  =  0.219086362515982;
            nw(9)  =  0.149451349150581;
            nw(10) =  0.066671344308688;
        case 12
            np(1)  = -0.981560634246719;
            np(2)  = -0.904117256370475;
            np(3)  = -0.769902674194305;
            np(4)  = -0.587317954286617;
            np(5)  = -0.367831498998180;
            np(6)  = -0.125233408511469;
            np(7)  =  0.125233408511469;
            np(8)  =  0.367831498998180;
            np(9)  =  0.587317954286617;
            np(10) =  0.769902674194305;
            np(11) =  0.904117256370475;
            np(12) =  0.981560634246719;

            nw(1)  =  0.047175336386512;
            nw(2)  =  0.106939325995318;
            nw(3)  =  0.160078328543346;
            nw(4)  =  0.203167426723066;
            nw(5)  =  0.233492536538355;
            nw(6)  =  0.249147045813403;
            nw(7)  =  0.249147045813403;
            nw(8)  =  0.233492536538355;
            nw(9)  =  0.203167426723066;
            nw(10) =  0.160078328543346;
            nw(11) =  0.106939325995318;
            nw(12) =  0.047175336386512;
    end
    
    % Adjust quadrature points based on dimension
    Q = NaN(npt*npt,2);
    W = NaN(npt*npt,1);
    n = 1;
    for i = 1:npt
        for j = 1:npt
            Q(n,:) = [np(i) np(j)];
            W(n)   = nw(i)*nw(j);
            n = n+1;
        end
    end
    
elseif strcmp(type,'TRI') == 1
    switch npt
        case 1
            gp(1,:) = [0.333333333333333 0.333333333333333];
            
            gw(1,:) = 1.000000000000000;
        case 3
            gp(1,:) = [0.166666666666667 0.166666666666667];
            gp(2,:) = [0.666666666666667 0.166666666666667];
            gp(3,:) = [0.166666666666667 0.666666666666667];

            gw(1,:) = 0.333333333333333;
            gw(2,:) = 0.333333333333333;
            gw(3,:) = 0.333333333333333;
        case 4
            gp(1,:) = [0.200000000000000 0.200000000000000];
            gp(2,:) = [0.600000000000000 0.200000000000000];
            gp(3,:) = [0.333333333333333 0.333333333333333];
            gp(4,:) = [0.200000000000000 0.600000000000000];

            gw(1,:) =  0.520833333333333;
            gw(2,:) =  0.520833333333333;
            gw(3,:) = -0.562500000000000;
            gw(4,:) =  0.520833333333333;
        case 6
            gp(1,:) = [0.091576213509771 0.091576213509771];
            gp(2,:) = [0.445948490915965 0.108103018168070];
            gp(3,:) = [0.816847572980459 0.091576213509771];
            gp(4,:) = [0.108103018168070 0.445948490915965];
            gp(5,:) = [0.445948490915965 0.445948490915965];
            gp(6,:) = [0.091576213509771 0.816847572980459];

            gw(1,:) = 0.109951743655322;
            gw(2,:) = 0.223381589678011;
            gw(3,:) = 0.109951743655322;
            gw(4,:) = 0.223381589678011;
            gw(5,:) = 0.223381589678011;
            gw(6,:) = 0.109951743655322;
        case 7
            gp(1,:) = [0.101286507323456 0.101286507323456];
            gp(2,:) = [0.470142064105115 0.059715871789770];
            gp(3,:) = [0.797426985353087 0.101286507323456];
            gp(4,:) = [0.333333333333333 0.333333333333333];
            gp(5,:) = [0.059715871789770 0.470142064105115];
            gp(6,:) = [0.470142064105115 0.470142064105115];
            gp(7,:) = [0.101286507323456 0.797426985353087];

            gw(1,:) = 0.125939180544827;
            gw(2,:) = 0.132394152788506;
            gw(3,:) = 0.125939180544827;
            gw(4,:) = 0.225030000300000;
            gw(5,:) = 0.132394152788506;
            gw(6,:) = 0.132394152788506;
            gw(7,:) = 0.125939180544827;
        case 9
            gp(1,:) = [0.165409927389841 0.037477420750088];
            gp(2,:) = [0.797112651860071 0.037477420750088];
            gp(3,:) = [0.437525248383384 0.124949503233232];
            gp(4,:) = [0.037477420750088 0.165409927389841];
            gp(5,:) = [0.797112651860071 0.165409927389841];
            gp(6,:) = [0.124949503233232 0.437527248383384];
            gp(7,:) = [0.437525248383384 0.437525248383384];
            gp(8,:) = [0.037477420750088 0.797112651860071];
            gp(9,:) = [0.165409927389841 0.797112651860071];

            gw(1,:) = 0.063691414286223;
            gw(2,:) = 0.063691414286223;
            gw(3,:) = 0.205950504760887;
            gw(4,:) = 0.063691414286223;
            gw(5,:) = 0.063691414286223;
            gw(6,:) = 0.205950504760887;
            gw(7,:) = 0.205950504760887;
            gw(8,:) = 0.063691414286223;
            gw(9,:) = 0.063691414286223;
        case 12
            gp(1,:)  = [0.063089014491502 0.063089014491502];
            gp(2,:)  = [0.310352451033785 0.053145049844816];
            gp(3,:)  = [0.636502499121399 0.053145049844816];
            gp(4,:)  = [0.873821971016996 0.063089014491502];
            gp(5,:)  = [0.249286745170910 0.249286745170910];
            gp(6,:)  = [0.501426509658179 0.249286745170910];
            gp(7,:)  = [0.053145049844816 0.310352451033785];
            gp(8,:)  = [0.636502499121399 0.310352451033785];
            gp(9,:)  = [0.249286745170910 0.501426509658179];
            gp(10,:) = [0.053145049844816 0.636502499121399];
            gp(11,:) = [0.310352451033785 0.636502499121399];
            gp(12,:) = [0.063089014491502 0.873821971016996];

            gw(1,:)  = 0.050844906370207;
            gw(2,:)  = 0.082851075618374;
            gw(3,:)  = 0.082851075618374;
            gw(4,:)  = 0.050844906370207;
            gw(5,:)  = 0.116786275726379;
            gw(6,:)  = 0.116786275726379;
            gw(7,:)  = 0.082851075618374;
            gw(8,:)  = 0.082851075618374;
            gw(9,:)  = 0.116786275726379;
            gw(10,:) = 0.082851075618374;
            gw(11,:) = 0.082851075618374;
            gw(12,:) = 0.050844906370207;
        case 13
            gp(1,:)  = [0.065130102902216 0.065130102902216];
            gp(2,:)  = [0.312865496004875 0.048690315425316];
            gp(3,:)  = [0.638444188569809 0.048690315425316];
            gp(4,:)  = [0.869739794195568 0.065130102902216];
            gp(5,:)  = [0.048690315425316 0.312865496004875];
            gp(6,:)  = [0.260345966079038 0.260345966079038];
            gp(7,:)  = [0.479308067841923 0.260345966079038];
            gp(8,:)  = [0.638444188569809 0.312865496004875];
            gp(9,:)  = [0.333333333333333 0.333333333333333];
            gp(10,:) = [0.260345966079038 0.479308067841923];
            gp(11,:) = [0.048690315425316 0.638444188569809];
            gp(12,:) = [0.312865496004875 0.638444188569809];
            gp(13,:) = [0.065130102902216 0.869739794195568];

            gw(1,:)  =  0.053347235608839;
            gw(2,:)  =  0.077113760890257;
            gw(3,:)  =  0.077113760890257;
            gw(4,:)  =  0.053347235608839;
            gw(5,:)  =  0.077113760890257;
            gw(6,:)  =  0.175615257433204;
            gw(7,:)  =  0.175615257433204;
            gw(8,:)  =  0.077113760890257;
            gw(9,:)  = -0.149570044467670;
            gw(10,:) =  0.175615257433204;
            gw(11,:) =  0.077113760890257;
            gw(12,:) =  0.077113760890257;
            gw(13,:) =  0.053347235608839;
    end
    
    Q = gp;
    W = gw;
end