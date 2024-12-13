function maskCell = makeCellMask (xImSz, yImSz, rCell,lCell)

xCellCent = ceil(xImSz/2)
yCellCent = ceil(yImSz/2)

cl = lCell-2*rCell;
ytop = yCellCent + rCell;
ybot =yCellCent - rCell;
xtop=xCellCent+cl/2+rCell;
xbot=xCellCent-cl/2-rCell;



maskCell = zeros(yImSz,xImSz);
maskCell (ybot:ytop,xCellCent-cl/2:xCellCent+cl/2)= 1

circle =zeros(floor(rCell)*2+1,floor(rCell)*2+1)'

[coordsCircX,coordsCircY] = meshgrid(-rCell:rCell, -rCell:rCell);
rPos = sqrt(coordsCircX.^2+coordsCircY.^2);

for iiX = 1:size(rPos,1)
    for iiY = 1:size(rPos,2)
        if rPos(iiX,iiY)<rCell
            circle(iiX,iiY)=1;
        end
    end
end

circR = zeros(size(circle))
circL = circR;
circR(:,rCell+1:end)=circle(:,rCell+1:end)
circL(:,1:rCell+1)=circle(:,1:rCell+1)

maskCell(ybot:ytop,xbot:xCellCent-cl/2+rCell)=maskCell(ybot:ytop,xbot:xCellCent-cl/2+rCell)|circL;
maskCell(ybot:ytop,xCellCent+cl/2-rCell:xtop)=maskCell(ybot:ytop,xCellCent+cl/2-rCell:xtop)|circR;
%figure; imshow(maskCell,[])

end
