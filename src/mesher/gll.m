function [ksi, wi] = gll(n)
    % 
    % Gauss-Lobatto-Legendre (GLL) points and corresponding weights
    % used for local coordinates of nodes for 1D spectral element
    % Note that local coordinates and weights are hard coded based on respective equations
    % 
    % USAGE:: 
    %
    %   [ksi, wi] = gll(n)
    % 
    % Arguments: 
    %     n (integer): 
    %       number of nodes in 1D element, element approximation order is (n-1)
    % 
    % Returns: 
    %     ksi (double): 
    %       local coordinates in range -1:1, row vector of dimensions [1,n]
    %
    %     wi (double): 
    %       weights, row vector of dimensions [1,n]
    %
    % 
    % (C) Copyright 2024 Pawel Kudela, pk@imp.gda.pl 
    % Institute of Fluid Flow Machinery, Polish Academy of Sciences 
    % Mechanics of Intelligent Structures Department 

%--------------------------------------------------------------------------------------------------- 

ksi = zeros(1, n);
wi = zeros(1, n);
switch n
    case 10
        ksi(1) = -1.;
        ksi(2) = -0.9195339081664589;
        ksi(3) = -0.738773865105505;
        ksi(4) = -0.4779249498104445;
        ksi(5) = -0.16527895766638692;
        ksi(6) = -ksi(5);
        ksi(7) = -ksi(4);
        ksi(8) = -ksi(3);
        ksi(9) = -ksi(2);
        ksi(10) = 1.;
        wi(1) = 0.022222222222222223;
        wi(2) = 0.1333059908510802;
        wi(3) = 0.22488934206312652;
        wi(4) = 0.29204268367968383;
        wi(5) = 0.32753976118389744;
        wi(6) = wi(5);
        wi(7) = wi(4);
        wi(8) = wi(3);
        wi(9) = wi(2);
        wi(10) = wi(1);
    case 9
        ksi(1) = -1.;
        ksi(2) = -0.8997579954114602;
        ksi(3) = -0.6771862795107377;
        ksi(4) = -0.36311746382617816;
        ksi(5) = 0.;
        ksi(6) = -ksi(4);
        ksi(7) = -ksi(3);
        ksi(8) = -ksi(2);
        ksi(9) = 1.;
        wi(1) = 0.027777777777777776;
        wi(2) = 0.16549536156080175;
        wi(3) = 0.27453871250016254;
        wi(4) = 0.34642851097304594;
        wi(5) = 0.37151927437641724;
        wi(6) = wi(4);
        wi(7) = wi(3);
        wi(8) = wi(2);
        wi(9) = wi(1);
    case 8
        ksi(1) = -1.;
        ksi(2) = -0.8717401485096066;
        ksi(3) = -0.5917001814331423;
        ksi(4) = -0.2092992179024788;
        ksi(5) = -ksi(4);
        ksi(6) = -ksi(3);
        ksi(7) = -ksi(2);
        ksi(8) = 1.;
        wi(1) = 0.03571428571428571;
        wi(2) = 0.2107042271435067;
        wi(3) = 0.3411226924835044;
        wi(4) = 0.4124587946587037;
        wi(5) = wi(4);
        wi(6) = wi(3);
        wi(7) = wi(2);
        wi(8) = wi(1);
    case 7
        ksi(1) = -1.;
        ksi(2) = -0.8302238962785669;
        ksi(3) = -0.4688487934707142;
        ksi(4) = 0.;
        ksi(5) = -ksi(3);
        ksi(6) = -ksi(2);
        ksi(7) = 1.;
        wi(1) = 0.047619047619047616;
        wi(2) = 0.2768260473615661;
        wi(3) = 0.43174538120986283;
        wi(4) = 0.4876190476190476;
        wi(5) = wi(3);
        wi(6) = wi(2);
        wi(7) = wi(1);
    case 6
        ksi(1) = -1.;
        ksi(2) = -0.7650553239294647;
        ksi(3) = -0.28523151648064504;
        ksi(4) = -ksi(3);
        ksi(5) = -ksi(2);
        ksi(6) = 1.;
        wi(1) = 0.06666666666666667;
        wi(2) = 0.37847495629784766;
        wi(3) = 0.5548583770354865;
        wi(4) = wi(3);
        wi(5) = wi(2);
        wi(6) = wi(1);
    case 5
        ksi(1) = -1.;
        ksi(2) = -0.6546536707079771;
        ksi(3) = 0.;
        ksi(4) = -ksi(2);
        ksi(5) = 1.;
        wi(1) = 0.1;
        wi(2) = 0.5444444444444445;
        wi(3) = 0.7111111111111111;
        wi(4) = wi(2);
        wi(5) = wi(1);
    case 4
        ksi(1) = -1.;
        ksi(2) = -0.4472135954999579;
        ksi(3) = -ksi(2);
        ksi(4) = 1.;
        wi(1) = 0.16666666666666666;
        wi(2) = 0.8333333333333333;
        wi(3) = wi(2);
        wi(4) = wi(1);
    case 3
        ksi(1) = -1.;
        ksi(2) = 0.0;
        ksi(3) = 1.;
        wi(1) = 0.3333333333333333;
        wi(2) = 1.3333333333333333;
        wi(3) = wi(1);
    case 2
        ksi(1) = -1;
        ksi(2) = 1;
        wi(1) = 1;
        wi(2) = 1;
    case 1
        ksi(1) = 0;
        w(1) = 1;
    otherwise 
        disp('Chosen number of nodes is not available'); 
        w=[]; 
        ksi=[];
        return;
end
%---------------------------------------------------------------------------------------------------

end


