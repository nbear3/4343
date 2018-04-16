function airfoils = load_airfoils()
    % load airfoils from the \airfoils directory
    files = dir('airfoils/*.xlsx');
    r = size(files);
    airfoils = cell(1, r);
    for i = 1:r/2
        airfoils(i).root = xlsread(strcat('airfoils/', files(2*i-1).name));
        airfoils(i).tip = xlsread(strcat('airfoils/', files(2*i).name));
    end
