rho0 = 1000;
Nx = 1286;
Ny = 936;
Nz = 128;
saveDir = fullfile(pwd,"random_media_new");
mkdir(saveDir)

rng('shuffle')
for aa = 1:8
    disp(aa)
    density = rho0*(1 + 0.02*randn(Nx,Ny,Nz, "single"));
    save(fullfile(saveDir,"BAPW_STD2_REF2025_"+aa+".mat"), "density")
end


rng('shuffle')
for aa = 1:4
    disp(aa)
    density = rho0*(1 + 0.02*randn(Nx,Ny,Nz, "single"));
    save(fullfile(saveDir,"BAPW_STD2_SAM2025_"+aa+".mat"), "density")
end