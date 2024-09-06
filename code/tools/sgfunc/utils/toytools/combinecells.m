function Combinedcell = combinecells(Cell1, Cell2)
%Combinedcell = combinecells(Cell1, Cell2)
%   returns a cell column-vector, of which each cell contains a horizontal concatenation of elements of cell1 and cell2 

N1 = numel(Cell1);
N2 = numel(Cell2);
Combinedcell = cell(N1*N2,1);
for i1 = 1:N1
    for i2 = 1:N2
        j = (i1-1)*N2 + i2;
        Combinedcell{j,1} = [Cell1{i1},Cell2{i2}];
    end
end

end