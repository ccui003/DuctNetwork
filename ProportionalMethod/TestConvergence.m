duct = DuctNetwork();
duct.AddParameter({'Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'},[1.204,9e-5,15.11e-6]);
duct.AddBranch(0,1);
duct.AddFitting(1,'FanQuadratic',[100,10]);%fan
duct.AddBranch(1,0);
duct.AddFitting(2,'CircularDarcyWeisbach',[4.6,0.3]);%duct1

duct.Sim();
