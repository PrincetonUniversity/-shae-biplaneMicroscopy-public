% Create a composite image from the two halves of a biplane image. 
function [compositeImage,newImage] = biplaneCompImage(filteredImage, transf, limits, imageSize)
    newImage= zeros(512,512,2);
    newImage(limits(2,2)+1:imageSize ,limits(1,1):limits(3,1),1)= filteredImage (limits(2,2)+1:imageSize ,limits(1,1):limits(3,1));
    for iiRow=1:limits(2,2)
        for iiCol=1:limits(3,1)
        avals=fnval(transf,[iiCol,iiRow]');
        newCol = round(avals(1));
        newRow = round(avals(2));
            if (newRow >0 && newCol >0 )
                newImage(newRow,newCol,2)=filteredImage(iiRow,iiCol);
            end
        end
    end
    compositeImage = newImage(:,:,1)+newImage(:,:,2);
    max(max(compositeImage));
    min(min(compositeImage)) ;
    %compositeImage = newImage(:,:,1)*newImage(:,:,2);