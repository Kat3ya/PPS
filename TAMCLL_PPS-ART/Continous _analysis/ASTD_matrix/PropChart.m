function Cprop=PropChart(Basins_size,Mc_all_sort)
C=mat2cell(Mc_all_sort,Basins_size,Basins_size');
for i1=1:size(C,1)
    for i2=1:size(C,2)
        %Cprop(i1,i2)=(1/length(Mc_all_sort)^2)*nnz(C{i1,i2})*100;
        Cprop(i1,i2)=(nnz(C{i2,i1})/nnz(Mc_all_sort))*100;
    end
end