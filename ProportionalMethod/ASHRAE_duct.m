duct = DuctNetwork();
duct.AddParameter({'Density(kg/m^3)','Roughness(mm)','Kinematic Viscosity(m^2/s)'},[1.204,9e-5,15.11e-6]);
duct.AddBranch(0,1);
duct.AddFitting(1,'CircularDarcyWeisbach',[4.6,0.3]);%duct1
duct.AddFitting(1,'ED1_3',[0.3,0.06]);%1
duct.AddFitting(1,'CD9_1',[0.3,0.3,0]);%2 damper 1
duct.AddBranch(0,1);
duct.AddFitting(2,'CircularDarcyWeisbach',[18.3,0.2]);%duct2
duct.AddFitting(2,'ED1_1',[0.2,0.00161,0]);%4
duct.AddFitting(2,'CD6_1',[0.2,1,0.7]);%4 screen
duct.AddFitting(2,'CD3_6',[0.2]);%5
duct.AddFitting(2,'CD9_1',[0.2,0.2,0]);%6 damper 2
duct.AddBranch(1,3);
duct.AddFitting(3,'CircularDarcyWeisbach',[6.1,0.3]);%duct3
duct.AddFitting(3,'CD9_1',[0.3,0.3,0]);%7 
duct.AddFitting([1,3,2],'CircularY30Junction',[0.3,0.3,0.2]);%3 ED5_1
duct.AddBranch(0,2);
duct.AddFitting(4,'RectangularDarcyWeisbach',[1.5,0.6,0.6]);%duct4
duct.AddFitting(4,'CR9_4',[0.6,0.6]);%10
duct.AddFitting(4,'ER4_3',[0.6,0.6,0.38,0.75]);%11
duct.AddFitting(4,'Louver',-25);%9
duct.AddBranch(2,3);
duct.AddFitting(5,'CircularDarcyWeisbach',[18.3,0.38]);%duct5
duct.AddFitting(5,'CD3_17',[0.38]);%12
duct.AddFitting(5,'CD9_1',[0.38,0.38,0]);%13 damper 3
duct.AddBranch(3,4);
duct.AddFitting(6,'CircularDarcyWeisbach',[9.1,0.45]);%duct6
duct.AddFitting(6,'CD9_3',[0.45]);%14
duct.AddFitting(6,'CD3_9',[0.45]);%15
duct.AddFitting(6,'ED7_2',[0.45,1.5,0.9]);%ED7_2
duct.AddFitting([3,6,5],'CircularY45Junction',[0.3,0.45,0.38]);%8 ED5_2
duct.AddBranch(12,0);
duct.AddFitting(7,'RectangularDarcyWeisbach',[4.3,0.25,0.25]);%duct7
duct.AddFitting(7,'CR3_3',[0.25,0.25,0.7,90]);%16
duct.AddFitting(7,'CR9_1',[0.25,0.25,0]);%17 damper 4
duct.AddFitting(7,'Louver',-25);%43
duct.AddBranch(12,0);
duct.AddFitting(8,'RectangularDarcyWeisbach',[1.2,0.25,0.25]);%duct8
duct.AddFitting(8,'CR9_4',[0.25,0.25]);%18
duct.AddFitting(8,'Louver',-25);%44
duct.AddBranch(11,12);
duct.AddFitting(9,'RectangularDarcyWeisbach',[7.6,0.25,0.5]);%duct9
duct.AddFitting(9,'SR3_1',[0.25,0.5,0.4]);%20
duct.AddFitting([7,9,8],'RectangularTJunction',[0.25,0.25,0.25,0.5,0.25,0.25]);%19 SR5_13
duct.AddBranch(9,11);
duct.AddFitting(10,'RectangularDarcyWeisbach',[13.7,0.25,0.4]);%duct10
duct.AddFitting(10,'CR9_1',[0.25,0.4,0]);%21 damper 5
duct.AddFitting(10,'CR3_10',[0.25,0.4]);%22
duct.AddFitting(10,'CR3_6',[0.25,0.4,90]);%23
duct.AddBranch(10,0);
duct.AddFitting(11,'RectangularDarcyWeisbach',[3,0.25,0.25]);%duct11
duct.AddFitting(11,'CR9_1',[0.25,0.25,0]);%25 damper 6
duct.AddFitting(11,'SR2_1',[0.25,0.25]);%26
duct.AddBranch(10,0);
duct.AddFitting(12,'RectangularDarcyWeisbach',[6.7,0.25,0.25]);%duct12
duct.AddFitting(12,'CR9_1',[0.25,0.25,0]);%28 damper 7
duct.AddFitting(12,'SR2_5',[0.25,0.25,0.45,0.45,0.6]);%29
duct.AddBranch(9,10);
duct.AddFitting(13,'RectangularDarcyWeisbach',[10.7,0.25,0.35]);%duct13
duct.AddFitting(13,'CR9_1',[0.25,0.35,0]);%30
duct.AddFitting([11,12,13],'RectangularYJunction',[0.25,0.25,0.25,0.35]);%27 SR5_14
duct.AddBranch(7,9);
duct.AddFitting(14,'RectangularDarcyWeisbach',[4.6,0.25,0.66]);%duct14
duct.AddFitting(14,'CR9_1',[0.25,0.66,0]);%31
duct.AddFitting([14,13,10],'RectangularYJunction',[0.25,0.66,0.35,0.4]);%24 SR5_1
duct.AddBranch(8,0);
duct.AddFitting(15,'RectangularDarcyWeisbach',[12.2,0.15,0.2]);%duct15
duct.AddFitting(15,'CR3_1',[0.15,0.2,1.5,90]);%48
duct.AddFitting(15,'SR2_6',[0.15,0.2,0.5]);%33
duct.AddFitting(15,'CR9_1',[0.15,0.2,0]);%34 damper 8
duct.AddBranch(8,0);
duct.AddFitting(16,'RectangularDarcyWeisbach',[6.1,0.15,0.2]);%duct16
duct.AddFitting(16,'SR2_3',[0.2,0.15,20,0.56713]);%36
duct.AddFitting(16,'CR6_1',[0.15,0.2,2,0.8]);%36 screen
duct.AddFitting(16,'CR9_1',[0.15,0.2,0]);%37 damper 9
duct.AddBranch(7,8);
duct.AddFitting(17,'RectangularDarcyWeisbach',[4.2,0.15,0.25]);%duct17
duct.AddFitting(17,'CR9_1',[0.15,0.25,0]);%38
duct.AddFitting([17,15,16],'RectangularYJunction',[0.15,0.25,0.2,0.2]);%35 SR5_1
duct.AddBranch(6,7);
duct.AddFitting(18,'RectangularDarcyWeisbach',[7,0.25,0.8]);%duct18
duct.AddFitting(18,'CR6_4',[0.25,0.8,0.025,0]);%39
duct.AddFitting(18,'SR4_1',[0.25,0.8,0.45,0.45]);%40
duct.AddFitting(18,'CR3_17',[0.8,0.25,1]);%41
duct.AddFitting(18,'CR9_6',[0.25,0.8]);%45
duct.AddFitting([14,18,17],'RectangularTJunction',[0.25,0.66,0.25,0.8,0.15,0.25]);%32 SR5_13
duct.AddBranch(5,6);
duct.AddFitting(19,'RectangularDarcyWeisbach',[3.7,0.45,0.8]);%duct19
duct.AddFitting(19,'SR7_17',[0.45,0.3,0.45,0.8,1]);%42
duct.AddFitting(19,'CR9_4',[0.45,0.8]);%47
duct.AddFitting(19,'Louver',-15);%46
duct.AddBranch(4,5);
duct.AddFitting(20,'FanQuadratic',[1e8,1.9]);%fan

duct.d_Sidx = [10, 22, 46, 66, 91, 102, 110, 148, 162];
duct.gdr_Bidx = [1,2,4,7,8,11,12,15,16];