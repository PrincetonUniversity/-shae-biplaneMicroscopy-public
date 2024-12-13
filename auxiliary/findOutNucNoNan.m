function [outNuc,lenTrajs,inPole] =findOutNucNoNan(struct,flagSize)
%%

if nargin==1
    flagSize = 0;
end
%%
goodTracks= struct.goodTracks;
dataVars = struct.dataVars;
sizeF= struct.sizeF;
%figure
  clear outNuc
outNuc = nan(sizeF,length(goodTracks))  
  
for iiCell=1:length(goodTracks)%[1:24,26:length(goodTracks)] %counter
    %close all
    %%
    iiCell
     try
   
    load (char(dataVars(iiCell)),'saveFileGrainRot','jFile','imageStruct');
    load (char(dataVars(iiCell)),'minBot','maxBot','cellMask','th');
    load (char(dataVars(iiCell)),'centRot','measurements','bw','Inuc');

    
 
%     jFile
%     
%     %Finding original box orientation from user selected cells
%     a=imageStruct(jFile).cells;
%     x= cat(1,a(1:2,1,iiCell),a(3:4,1,iiCell));
%     y= cat(1,a(1:2,2,iiCell),a(3:4,2,iiCell));
%     bw=poly2mask(x,y,512,512);
%     figure; imshow(bw,[])
%     se90 = strel('line',2,90);
%     se0 = strel('line',2,0);
%     bw = imdilate(bw,[se90 se0]); %close edges
%     figure; imshow(bw,[])
% 
% 
%     I2 = uint8(255 * mat2gray(Inuc));
%     figure; imshow(double(I2).*double(bw),[])
% 
%     [~,thresh] =edge(I2, 'log');
%     edges =edge(I2, 'log',1*thresh); %seems to work best for nucleoids
%     edges = bw.*edges;
%     figure; imshow(edges,[]) 
% 
%     se90 = strel('line',1,90);
%     se0 = strel('line',1,0);
%     BWsdil4 = imdilate(edges,[se90 se0]); %close edges
%     figure; imshow(BWsdil4,[])
% 
%     BWsdil2 =filledgegaps(BWsdil4,3);
%     figure; imshow(BWsdil2,[])
%     BWdfill = imfill(BWsdil2,'holes'); %fill contours
%     figure; imshow(BWdfill,[])
% 
%     masked =BWdfill;%bw.*BWdfill;
% 
%     bw = bwareaopen(masked, 15);     %Removes noise 
%     figure;imshow(bw,[])
% 
%      BWoutline = bwperim(bw);
%         Segout = Inuc; 
%         Segout(BWoutline) = 255; 
%        figure, imshow(Segout,[],'initialMagnification','fit'), title('outlined original image');
%         %hold on
% 
%     measurements =regionprops(bw,'centroid','minorAxisLength','majorAxisLength','Extrema','Orientation','BoundingBox');
% 
% 
% 
%     figure; imshow(bw,[])
%     hold on
%     plot(measurements.Centroid(1),measurements.Centroid(2),'*')

    
  
    
    
    rotMask  =imrotate(cellMask,-th,'bicubic');   


    nucMask=bw(minBot(2):maxBot(2),minBot(1):maxBot(1));
    InucMask = Inuc(minBot(2):maxBot(2),minBot(1):maxBot(1));
    
    nucRot = imrotate(nucMask,-measurements.Orientation, 'bicubic');
    nucRotIm = imrotate(InucMask,-measurements.Orientation, 'bicubic');
   
     %% figure;imshow(double(nucRotIm),[])
     BWoutline = bwperim(nucMask);
    Segout = InucMask; 
    Segout(BWoutline) = 255; 
  % figure, imshow(Segout,[],'initialMagnification','fit'), title('outlined original image');
   % hold on
    
 
    %nucRot(bwperim(nucRot))=0;
   %  figure;imshow(nucRot)
   nucMeas=regionprops(nucRot,'centroid','minorAxisLength','majorAxisLength','Extrema','Orientation','BoundingBox');

cX=centRot(1);
cY=centRot(2);
cZ=25.75;
xDisp=nucMeas.Centroid(1)-cX;
yDisp=nucMeas.Centroid(2)-cY;
zDisp=nucMeas.Centroid(2)-cZ;

%%

%%
clear inXZ
clear inXY
nucTransXY =imtranslate(nucRot,[-xDisp,-yDisp],'FillValues',0);

%figure; imshow(nucTransXY,[])






%
bound= bwboundaries(nucTransXY);
boundy=bound{1};
yNuc=boundy(:,2);
xyNuc=boundy(:,1);
% hold on
% plot(boundy(:,2),boundy(:,1))

nucTransXZ =imtranslate(nucRot,[-xDisp,-zDisp],'FillValues',0);
%figure; imshow(nucTransXZ,[])
bound= bwboundaries(nucTransXZ);
boundy=bound{1};
zNuc=boundy(:,2);
xzNuc=boundy(:,1);
% hold on
% plot(boundy(:,2),boundy(:,1))





load (char(dataVars(iiCell)),'centRot','cropImTop','cropImBot');
load(char(dataVars(iiCell)),'Irot');
traj=csvread(goodTracks{iiCell});

%%
temp = traj(:,3);
for iiL= 1:length(traj(:,3))
    if isnan(traj(iiL, 3))
        if iiL ==1
            before = nan;
        else
            before = traj(iiL-1,3);
        end
        if iiL == length (traj(:,3))
            after= nan;
        else
            after = traj(iiL+1,3);
        end
        
        
        cMin =0;
        cMax =0;
        while isnan(before)
            before = traj(iiL-cMin,3);
            cMin =cMin+1;
            if (iiL-cMin==0)
                break
            end
        end
        while isnan(after)
            after = traj(iiL+cMax,3);
            cMax = cMax +1;
            if (iiL+cMax>=length(traj))
                break
            end
        end
        if iiL==1
            temp(iiL)=after;
        elseif iiL ==length(traj(:,3))
            temp(iiL)=before;
        else 
            temp(iiL)=(before+after)/2;
        end
    end
end
traj(:,3)=temp;
%%
lenTrajs(iiCell)=length(traj);

if (flagSize)
    sizeF=length(traj);
end
xTraj=traj(1:sizeF,1);
yTraj=traj(1:sizeF,2);
zTraj=(traj(1:sizeF,3)-cZ)/0.59+cZ;

     
[leftEnd,rightEnd] = slopeChangeEndsEdges(nucTransXY,1,1);

%leftEnd = min(zNuc);
%rightEnd = max(zNuc);

inPoleTemp = xTraj>=rightEnd | xTraj<=leftEnd;
%xTraj(inPole);
inPole(1:sizeF,iiCell)= inPoleTemp;
[inXZ,onXZ] = inpolygon(xTraj,zTraj,zNuc,xzNuc);
[inXY,onXY] = inpolygon(xTraj,yTraj,yNuc,xyNuc);

%%figure ; imshow(nucTransXZ)
 xv = xzNuc; yv = zNuc;
        xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
        x = xTraj; y = zTraj;
        in = inpolygon(x,y,xv,yv);
      %%  hold on; plot(zNuc,xzNuc,x(in),y(in),'.r',x(~in),y(~in),'.b')
%  


outNuc(1:sizeF,iiCell) = (~inXZ | ~inXY);

 test=(~inXZ | ~inXY)
% figure;plot(yNuc,xyNuc,'b-','linewidth',1.5)
% hold on 
% plot(traj(:,1),traj(:,2),'kx')

%plot(zNuc,xzNuc,'r-','linewidth',1.5)
%hold on 
%plot(traj(1:sizeF,1),(traj(1:sizeF,3)-cZ)/0.59+cZ,'kx')


%%
test =inPoleTemp;
% figure;plot(yNuc,xyNuc,'k-','linewidth',1.5)
% hold on 
% plot(traj(test,1),traj(test,2),'.k')
% plot(traj(~test,1),traj(~test,2),'m.','markersize',8)
% 
% 
% plot(zNuc,xzNuc,'r-','linewidth',1.5)
% hold on 
% plot(traj(1:sizeF,1),(traj(1:sizeF,3)-cZ)/0.59+cZ,'kx')



%%

   
   %% FIgure of one trajectory
   %[a,b]= uigetfile();
%    traj=csvread(    'C:\Users\dsmendez\Documents\GitHub\shae-miscTools\Diana\20200911_charge_exponential\cropped\+9\20200911_174216.75\5\5_traj_cell_1.csv');
%    
%    
%    
%      figure;plot3((traj(:,1))*80,(traj(:,2))*80,(traj(:,3))/0.59,'k.-')
%      set(gca, 'FontSize', 20)
%      axis equal

     catch
        warning('error')
     end
    
    
    
    

end
