

files = dir('airfoils/*.xlsx');
r = size(files);
airfoil = cell(1, r);
for i = 1:r/2
    airfoil(i).tip = xlsread(strcat('airfoils/', files(2*i-1).name));
    airfoil(i).root = xlsread(strcat('airfoils/', files(2*i).name));
end
