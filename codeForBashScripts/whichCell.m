function [jFile,iiCell]= whichCell(imageStructPath, t)
load(imageStructPath)
for iiCell =1:size(imageStruct,2)
   sizes(iiCell)= size(imageStruct(iiCell).cells,3)
end
sumSizes = cumsum(sizes);
limits =[1, cumsum(sizes)+1]
% % 
jFile = discretize(t,limits)


if( jFile ==1)
    iiCell = t
else
    iiCell = t- sumSizes(jFile-1)
end

end