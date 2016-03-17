function GBF = generateMagneticFieldFile(pathToFile,spacing,NP)
% Function to generate a table of values of the magnetic field components
% in cylindrical coordinates (2D).
% EXAMPLES: 
% generateMagneticFieldFile('ANALYTIC2D.dat','equidistant',[150,150])
% generateMagneticFieldFile('PADUA.dat','padua',150^2)
% generateMagneticFieldFile('CHEBYSHEV.dat','chebyshev',[50,50])


narginchk(3,3)

B = analyticalB([],'initialize');

offset = 2;

Rmin = (B.Ro - B.a) - offset;
Rmax = (B.Ro + B.a) + offset;
Zmin = - (B.a + offset);
Zmax = B.a + offset;

switch spacing
    case 'equidistant'
        NR = NP(1);
        NZ = NP(2);
        R = linspace(Rmin,Rmax,NR);
        Z = linspace(Zmin,Zmax,NZ);
    case 'padua'
        n = 1; % initial guess for order of padua points
        N = (n + 1)*(n + 2)/2;
        while (N - NP < 0)
            n = n + 1;
            N = (n + 1)*(n + 2)/2;
        end
        RZ = paduapts(n, [Rmin, Rmax, Zmin, Zmax]);
        R = RZ(:,1);
        Z = RZ(:,2);
        
        disp(['Using Padua points of degree ' num2str(n)])
        disp(['Array dimension [' num2str(N) ',' num2str(N) ']'])
    case 'chebyshev'
        NR = NP(1);
        NZ = NP(2);
        [R,Z] = chebpts2(NR,NZ,[Rmin, Rmax, Zmin, Zmax]);
    otherwise
        error('Enter a valid option.')
end


if strcmp(spacing,'equidistant')
    
elseif strcmp(spacing,'padua')
   
end

data = zeros(1,7);
fileID = fopen(pathToFile,'w');
switch spacing
    case 'equidistant'
        for ix=1:NR
            for iz=1:NZ
                [data(5:7),data(1:3)] = analyticalB([R(ix),0,Z(iz)]);
                fprintf(fileID,'%3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\n',data);
            end
        end
    case 'padua'
        for ii=1:N
            [data(5:7),data(1:3)] = analyticalB([R(ii),0,Z(ii)]);
            fprintf(fileID,'%3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\n',data);
        end
    case 'chebyshev'
        for ix=1:NR
            for iz=1:NZ
                [data(5:7),data(1:3)] = analyticalB([R(ix,iz),0,Z(ix,iz)]);
                fprintf(fileID,'%3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\n',data);
            end
        end
end

fclose(fileID);

disp(['File: ``' pathToFile '`` generated!'])
end

function [B,P] = analyticalB(X,opt)
% Analytical magnetic field
% X is a vector X(1)=R, X(2)=y, X(3)=Z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Bo = 5;
a = 2;% 0.6;% Minor radius in meters.
Ro = 6.0; % Major radius in meters.
qa = 10; % Safety factor at the separatrix (r=a)
co = 0.5; % Extra parameter
lamb = a/co;
Bpo = (a/Ro)*(Bo/qa)*(1-co^2)/co;
% Parameters of the analytical magnetic field

if nargin == 2
    if strcmp(opt,'initialize')
        B = struct;
        B.Bo = Bo;
        B.a = a;% 0.6;% Minor radius in meters.
        B.Ro = Ro; % Major radius in meters.
        B.qa = qa; % Safety factor at the separatrix (r=a)
        B.co = co; % Extra parameter
        B.lamb = lamb;
        B.Bpo = Bpo;
    elseif strcmp(opt,'normalize')
    end
else
    % Toroidal coordinates
    % r = radius, theta = poloidal angle, zeta = toroidal angle
    r = sqrt( (sqrt(X(1)^2 + X(2)^2) - Ro)^2 + X(3)^2 );
    theta = atan2(X(3),sqrt(X(1)^2 + X(2)^2) - Ro);
    if theta < 0
        theta = theta + 2*pi;
    end
    zeta = atan2(X(1),X(2));
    if zeta < 0
        zeta = zeta + 2*pi;
    end
    % Toroidal coordinates
    
    % Poloidal magnetic field
    Bp = Bpo*(r/lamb)/( 1 + (r/lamb)^2 );
    
    eta = r/Ro;
    Br = 1/( 1 + eta*cos(theta) );
    
    Bx = Br*( Bo*cos(zeta) - Bp*sin(theta)*sin(zeta) );
    By = -Br*( Bo*sin(zeta) + Bp*sin(theta)*cos(zeta) );
    Bz = Br*Bp*cos(theta);
    
    % Now we change to cylindrical coordinates
    % Cylindrical coordinates
    % R = radius, phi = azimuthal angle, Z = Z coordinate
    R = sqrt(X(1)^2 + X(2)^2);
    Z = X(3);
    phi = atan2(X(2),X(1));
    if phi < 0
        phi = phi + 2*pi;
    end
    
    % B = [B^R,B^phi,B^Z]
    B = [Bx*cos(phi) + By*sin(phi),...
        By*cos(phi) - Bx*sin(phi),...
        Bz];
    P = [R,phi,Z];
end

end
