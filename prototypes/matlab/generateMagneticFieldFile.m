function GBF = generateMagneticFieldFile(pathToFile,res)
narginchk(2,2)

Nx = res(1);
Nz = res(2);

B = analyticalB([],'initialize');

xmin = (B.Ro - B.a) - 1;
xmax = (B.Ro + B.a) + 1;
zmin = - (B.a + 1);
zmax = B.a + 1;
x = linspace(xmin,xmax,Nx);
z = linspace(zmin,zmax,Nz);

numVals = Nx*Nz;
data = zeros(1,7);

fileID = fopen(pathToFile,'w');

for ix=1:Nx
    for iz=1:Nz
        [data(5:7),data(1:3)] = analyticalB([x(ix),0,z(iz)]);
        fprintf(fileID,'%3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\n',data);
    end
end
fclose(fileID)

disp(['File: ``' pathToFile '`` generated!'])
end

function [B,P] = analyticalB(X,opt)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

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
    % R = radius, phi = azimuthal angle, Z = z coordinate
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
