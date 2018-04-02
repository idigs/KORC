iSpp.spp1.A = 18; % Argon
iSpp.spp1.nn = 0.002611*1E13*1E6; % neutral density in m^-3
iSpp.spp1.In = 15.7596117; % Ionization energy of neutral Argon.
iSpp.spp1.nz = [0.284,1.17]*1E13*1E6; % density of partially ionized impurities.
iSpp.spp1.Z = [1,2]; % Average charge state.
iSpp.spp1.Iz = [27.62967,40.735]; % Ionization energy.

iSpp.spp2.A = 1; % Deuterium
iSpp.spp2.nn = 0.04481*1E13*1E6; % neutral density in m^-3
iSpp.spp2.In = 15.46658; % Ionization energy of neutral Argon.
iSpp.spp2.nz = 3.764*1E13*1E6; % Ionized deuterium.
iSpp.spp2.Z = 1;
iSpp.spp2.Iz = [];