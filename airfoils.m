files = dir('airfoils/*.xlsx');
[r c] = size(files);
airfoil = {};
cd 'airfoils'
for i = 1:r
    name = files(i).name;
    % Odd is tip
    ind = floor(i/2) + mod(i,2);
    if mod(i, 2) == 1
        airfoil(ind).tip = xlsread(name);
    else
        % Even is root
        airfoil(ind).root = xlsread(name);
    end
end

cd ..
