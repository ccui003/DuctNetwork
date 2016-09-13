FileName='SR5_13.txt';
Cb2 = zeros(11,11);
Cb2 = Get2DTable(FileName, 6,16)';
QbQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AbAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
Cb = Cb2.*(QbQc'.^2*AbAc.^-2);

Cs2 = zeros(11,11);
Cs2 = Get2DTable(FileName, 23,33)';
QsQc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
AsAc = [0.01      0.1      0.2      0.3      0.4      0.5      0.6      0.7      0.8      0.9      1.0];
Cs = Cs2.*(QsQc'.^2*AsAc.^-2);

save('SR5_13.mat','Cb','Cb2','QbQc','AbAc','Cs','Cs2','QsQc','AsAc');