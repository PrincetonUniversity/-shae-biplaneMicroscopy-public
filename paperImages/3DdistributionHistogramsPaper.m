
%%
close all
index=1
sizeF=750
 clear resFull
for iiCell=[5,6,8,11,12,13,16]%1:length(goodTracks)
    %[1:18,20:length(goodTracks)]%
    %[1:10,12:length(goodTracks)]pfvfinal
   % [1,3,11,18] %splogged pfv
    %-18 [1,3:6,9:length(goodTracks)]
    %[1:18,20:length(goodTracks)] full +15
    %[1:6,8,10:13,15:17,20:length(goodTracks)]+15
   % [1:6,8,10:13,15:17,20:length(goodTracks)-1]
    %1:length(goodTracks) %aqls
    %[2,4:10,12:17,19:length(goodTracks)]%pfv
    %[1:6,7:9,10:15,17:length(goodTracks)-1] vuldi
    %[1:length(goodTracks)-1]
    %1:length(goodTracks)
    %[5,6,8,11,12,13,16]%1:length(goodTracks)
    %[2:6,8:10,12:15,17,19:length(goodTracks)]
    %[2,4:10,12:17,19:length(goodTracks)]
 %   index=1
%sizeF=500
 

%clear resFull
    
    %[1:6,8,9,11,12,13,15,16,19:length(goodTracks)]
    %-18 [1,3:6,9:15,17:length(goodTracks)]
    %pfv [2,4:10,12:length(goodTracks)]%
    %[1,3,5,12,13,17]%-18
    %8,12,
    iiCell  
    load (char(dataVars(iiCell)),'centRot','rotMatrix','Irot');
    load(char(dataVars(iiCell)),'cellMask','measurements');
    load (char(dataVars(iiCell)),'minBot','maxBot','Inuc');
    load (char(dataVars(iiCell)),'cropImTop','cropImBot','minBot');
    
    %cellMaskNot = ~cellMask
    %figure;imshow(cellMask,[])
       
    measurementsCell =regionprops(cellMask,'Orientation','BoundingBox','minorAxisLength','majorAxisLength');
    
     
    InucTest= Inuc(minBot(2):maxBot(2),minBot(1):maxBot(1));
    IMaskedTest= double(cellMask).*double(InucTest(:,:));
    imageIn =cropImTop(:,:,1)+ cropImBot(:,:,1);

    a0(1) = round(measurements.MajorAxisLength)+4; %3um
    a0(2) = round(measurements.MinorAxisLength/2)+1; %1um 5
    a0(3) = measurements.Centroid(1)-minBot(1)+1;
    a0(4) = measurements.Centroid(2)-minBot(2)+1;
    a0(5) = 0.7;
    a0(6) = measurements.Orientation;

   
    
    
    [fittedCellNoConv,cellParams,cellMask1] = cellFit2StageNoConvRotTest(imageIn,a0);

    
    
    BWoutline = bwperim(cellMask);
    Segout = InucTest; 
    Segout(BWoutline) = 255; 
    %figure, imshow(Segout,[],'initialMagnification','fit'), title('outlined original image');
    %hold on
    %rectangle('Pos',measurementsCell.BoundingBox)
    
    %lcell = paramsConv(1)+2+2*paramsConv(2);
    lcell = cellParams(1)-4;
    nucToCell(iiCell)=measurements.MajorAxisLength/lcell;
    dcell = cellParams(2)*2;
    dcellTolcell(iiCell)=dcell/lcell;
    traj=csvread(goodTracks{iiCell});
    %figure;plot(traj(:,2)*80,traj(:,3)*100/0.59);axis equal
    
    
    %drawnow
    sizeF = length(traj)
    centX= measurements.Centroid(1)-minBot(1)+1;
    centY = measurements.Centroid(2)-minBot(2)+1;
    
    centZero = [centX-size(InucTest,2)/2,centY-size(InucTest,1)/2];
    centRot1 =centZero * rotMatrix;
    centRot1 = [centRot1(:,1)+size(InucTest,2)/2 , centRot1(:,2)+size(InucTest,1)/2];
    centRot1 = [centRot1(:,1)+ (size(Irot,2)-size(InucTest,2))/2 , centRot1(:,2)+ (size(Irot,1)-size(InucTest,1))/2];
   
    
    x=(traj(1:sizeF,1)-centRot1(1))*30/lcell;

    y=(traj(1:sizeF,2)-centRot1(2))*10/dcell;

    centZ =  25.75%28;%-6% z0(iiCell)-2;25.75%
    %For testing circ trajs
    %%centZ = (max(traj(1:sizeF,3))-min(traj(1:sizeF,3)))/2+min(traj(1:sizeF,3));
    %centY = (max(traj(1:sizeF,2))-min(traj(1:sizeF,2)))/2+min(traj(1:sizeF,2));
    z = (traj(1:sizeF,3)-centZ)*10/dcell;
    figure;plot((traj(:,2)-centY)*80,(traj(:,3)-centZ)*100/0.59);axis equal
 
    resCent = [x,y,z,traj(1:sizeF,4:end)];
    
    
    resFull(index:index+length(resCent)-1,:) = resCent(:,1:6) ;   
    index = index +length(resCent);

end

%% After trajectories are all aligned 3D histograms are found and plotted.
% The individual panels were changed to black and white color scale and
% then the images were arranged as in figures 2.a, 3.a and S __


widthYZ=600
widthX=1800
xRes=resFull(:,1)*80;
yRes=-(resFull(:,2))*80;

%zRes=((resFull(:,3))+0.75)*100/0.59;%Aqls and -18
zRes=((resFull(:,3)))*100/0.59;

pixelSize =40

xyRes =[xRes50exp,yRes50exp]
zyRes=[yRes50exp,zRes50exp]
xzRes=[xRes50exp,zRes50exp]

edges ={(-widthX:pixelSize:widthX)',(-widthYZ:pixelSize:widthYZ)'};
figure
hist3(xyRes,'CDataMode','auto','FaceColor','interp','edges',edges)
xlabel('x')
ylabel('y')
title (strcat('cell',num2str(iiCell)))
axis equal
xlim([-widthX,widthX])
ylim([-widthYZ,widthYZ])
xticks([-widthX,-widthYZ,widthYZ,widthX])
yticks([-widthYZ,0,widthYZ])
set(gca, 'FontSize', 20)


edges ={(-widthYZ:pixelSize:widthYZ)',(-widthYZ:pixelSize:widthYZ)'}
figure
hist3(zyRes,'CDataMode','auto','FaceColor','interp','edges',edges)
xlabel('y')
ylabel('z')
title (strcat('cell',num2str(iiCell)))
axis equal
xlim([-widthYZ,widthYZ])
ylim([-widthYZ,widthYZ])
xticks([-widthYZ,0,widthYZ])
yticks([-widthYZ,0,widthYZ])
set(gca, 'FontSize', 20)


edges ={(-widthX:pixelSize:widthX)',(-widthYZ:pixelSize:widthYZ)'}
figure
hist3(xzRes,'CDataMode','auto','FaceColor','interp','edges',edges)
xlabel('x')
ylabel('Z')
title (strcat('cell',num2str(iiCell)))
axis equal
xlim([-widthX,widthX])
ylim([-widthYZ,widthYZ])
xticks([-widthYZ,0,widthYZ])
yticks([-widthYZ,0,widthYZ])
set(gca, 'FontSize', 20)




