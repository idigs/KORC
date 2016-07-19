function IT = calculate_params(B,Ro,a,Eo,eta,DT,numTransits)
% + B is the magnetic field
% + Ro is the major radius of the plasma
% + Eo is the range of energies in eV
% + eta is the pitch angle
% + DT is the time step as a fraction of the gyro-period
% + numbTransits is the number of transits around the torus of major radius
% Ro

kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
e0 = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass

eta = eta*pi/180;

gamma = (qe*Eo + me*c^2)/(me*c^2);
v = c*sqrt(1 - 1./gamma.^2);
vpar = v*cos(eta);

wc = qe*B./(gamma*me);
Tc = 2*pi./wc;

dt = DT*(2*pi./wc);

Ba = B*Ro/(Ro + a);
wca = qe*Ba./(gamma*me);
Tca = 2*pi./wca;

P = 2*pi*Ro;
simTime = numTransits*P./vpar;

IT = simTime./dt;
% IT = ceil(Tca./dt);
disp(['Eo = ' num2str(Eo/1E6) ' IT = ' num2str(IT)])

end