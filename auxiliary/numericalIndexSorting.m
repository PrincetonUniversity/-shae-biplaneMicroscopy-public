function sortedIdx = numericalIndexSorting(filenames, pattern )
clear order1
%pattern='cell_\d*'
order=regexp(filenames, pattern, 'match')
order=string(order)
order=regexp(order, '\d*', 'match')
for i=1:length(filenames)
    order1(i)=str2num(char(order{i}));
end
[~,sortedIdx]=sort(order1);