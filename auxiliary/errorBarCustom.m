function [d0,d1,d2] =errorBarCustom (xIn, yIn, yMin, yMax, lweb)
clear d1; clear d2

d0(1)=plot(xIn,yIn)
counter=2;
    for iiDot=1:length(xIn)
      
      %d2(iiDot)= plot([xIn(iiDot),xIn(iiDot)],[yMin(iiDot),yMax(iiDot)]);
      d1(iiDot)= plot([xIn(iiDot),xIn(iiDot)],[yMin(iiDot),yMax(iiDot)]);

      d2(counter*(iiDot-1)+1)= plot([xIn(iiDot)-lweb,xIn(iiDot)+lweb],[yMin(iiDot),yMin(iiDot)]);
      d2(counter*(iiDot-1)+2)= plot([xIn(iiDot)-lweb,xIn(iiDot)+lweb],[yMax(iiDot),yMax(iiDot)]);

    end

end