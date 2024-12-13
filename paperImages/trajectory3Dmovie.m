 % requires PATCHLINE from fileexchange or shae-miscTools
%positions =  csvread ('20190517_tridyes\msfGFP\20190517_185116.17_rotatedTraj.csv');
%positions =  csvread ('20201103_cccp_Everywhere\Cropped\pfvInduced\20201104_001805.83\5\5_traj_cell_2.csv');
%positions =  csvread ('20201103_cccp_Everywhere\Cropped\pfvInduced\20201103_233842.40\5\5_traj_cell_7.csv');
%positions =  csvread ('20201103_cccp_Everywhere\Cropped\pfvInduced\20201103_234000.24\5\5_traj_cell_7.csv');
%positions =csvread('20200730_32vs37\37at37\cropped\20200730_015656.34\5\5_traj_cell_2.csv')
%positions =  csvread ('20201103_cccp_Everywhere\Cropped\pfvInduced\20201103_234608.57\5\5_traj_cell_1.csv');
%positions =csvread('20191213_labeledVsUnlabeledOvernights\glucoseLabeled\20191213_230431.48\5\5_traj_cell_1.csv')
%positions = csvread('C:\Users\dsmendez\Documents\GitHub\shae-miscTools\Diana\20200911_charge_exponential\cropped\+9\20200911_174040.48\5\5_traj_cell_1.csv');
%positions= csvread ('C:\Users\dsmendez\Documents\GitHub\shae-miscTools\Diana\20200911_charge_exponential\cropped\+9\20200911_173414.86\5\5_traj_cell_7.csv');
positions = csvread (    'Z:\dsmendez\diana\code\Diana\20211014_ovChargeSweep\PFV\211014_204722timeseries\traj_cell_5_fixedZFull3.csv')
positions = csvread (    'Z:\dsmendez\diana\code\Diana\20211013_newFilterOvnights\AqLS\100msExp\211013_185640timeseries\traj_cell_18_fixedZFull3.csv')
positions = csvread(     'Z:\dsmendez\diana\code\Diana\20211014_ovChargeSweep\Vuldi\211015_013052timeseries\traj_cell_11_fixedZFull.csv')

%positions = csvread(    'Z:\dsmendez\diana\code\Diana\20211014_ovChargeSweep\Vuldi\211015_012314timeseries\traj_cell_1_fixedZFull.csv')

trajXVector =positions(:,1) *80;
trajYVector = positions(:,2)* 80;
trajZVector = positions(:,3) *100 ;
zeroX = centX*80%(max(trajXVector)+min(trajXVector))/2;
zeroY = centY*80%(max(trajYVector)+min(trajYVector))/2;
zeroZ = centZ*100%(max(trajZVector)+min(trajZVector))/2;
Xprime = trajXVector - zeroX;
Yprime = trajYVector - zeroY;
Zprime = (trajZVector - zeroZ)/0.59;
t = positions(:,5);
%%
 
% generate a diffusive trajectory
% t = 1:500;
% x = cumsum(randn(1,500));
% y = cumsum(randn(1,500));
% z = cumsum(randn(1,500));

%t= timeStamp;
x =Xprime';
y=Yprime';
z=Zprime';
%clf;
%scatter3(x,y,z);
%figure(gcf);
%axis vis3d;
%%
colorX = [55,126,184]/255; %yz
colorY =[77,175,74]/255;%xz
colorZ =[228,26,28]/255;%xy

%timeLagMax = 50; % prepend trajectory with this many nans, 
timeLagMax = 100
% as the plot object will expect this many coordinates;
xNanPadded = padarray(x,[0,timeLagMax],nan,'pre');
yNanPadded = padarray(y,[0,timeLagMax],nan,'pre');
zNanPadded = padarray(z,[0,timeLagMax],nan,'pre');

 
figure('Position',[300, 241 , 888 , 787],'Color','w');
hold on

box off
set(gca,'Position', [0.2,0.2,0.6,0.6])
axis equal;
xlim([-1800,1300]);
ylim([-800,800]);
zlim([-800,800]);

set(gca,'FontSize',24)
set(gca,'LineWidth',1.5)

view(34,20);
lighting('gouraud')
camlight('local')
camlight('left','local')
% view(90,0)
% camlight('headlight','local')
% view(0,90)
% camlight('headlight','infinite')

patchline(xlim',[min(ylim),min(ylim)]',[max(zlim),max(zlim)]','edgecolor','k','linewidth',1.5);
patchline(xlim,[max(ylim),max(ylim)],[max(zlim),max(zlim)],'edgecolor','k','linewidth',1.5);
patchline(xlim,[max(ylim),max(ylim)],[min(zlim),min(zlim)],'edgecolor',[0.7,0.7,0.7],'linewidth',1.5,...
    'linestyle','--','linewidth',1.5);

patchline([min(xlim),min(xlim)],ylim,[max(zlim),max(zlim)],'edgecolor','k','linewidth',1.5);
patchline([min(xlim),min(xlim)],ylim,[min(zlim),min(zlim)],'edgecolor',[0.7,0.7,0.7],...
    'linestyle','--','linewidth',1.5);
patchline([max(xlim),max(xlim)],ylim,[max(zlim),max(zlim)],'edgecolor','k','linewidth',1.5);


patchline([max(xlim),max(xlim)],[min(ylim),min(ylim)],zlim,'edgecolor','k','linewidth',1.5);
patchline([max(xlim),max(xlim)],[max(ylim),max(ylim)],zlim,'edgecolor','k','linewidth',1.5);
patchline([min(xlim),min(xlim)],[max(ylim),max(ylim)],zlim,'edgecolor',[0.7,0.7,0.7],...
    'linestyle','--','linewidth',1.5);

% curve = [xlim;[min(ylim),min(ylim)];[max(zlim),max(zlim)]];
% r=5;
% n=20;
% 
% tubeplot(curve,r,n,r/2,0)

ambStr=0.5;
diffStr= ambStr;
specStr = 1;
%%
%YZ (blue)
nX = 50; 
rCirc = 500 %sqrt(max(z.^2+y.^2));
theta = linspace (0,2*pi(), nX);
xCirc =min(xlim);
yCirc = rCirc*cos(theta)';
zCirc = rCirc*sin(theta)';

%xCircV=[xCirc;min(xlim)];
yCircV=[yCirc;0];
zCircV=[zCirc;0];
vCirc = [yCircV,yCircV, zCircV];
vCirc(:,1) = xCirc;

fCirc =cat(2,(1:nX-1)',(2:nX)', (nX+1)*ones(nX-1,1));
fCirc = [fCirc;[nX,1,nX+1]];

patch('Faces',fCirc,'Vertices',vCirc,'FaceColor',[0.85,0.85,0.85],'EdgeColor','none'...);
    ,'specularStrength',specStr,'diffuseStrength',diffStr,'ambientStrength',ambStr,'FaceAlpha',0.7);


%XZ (green)
xCylXZ = sign(yCircV)*500-double(yCircV<0)*500+yCircV;
yCylXZ = max(ylim);
zCylXZ = zCircV;
vCylXZ = [xCylXZ,xCylXZ, zCylXZ];
vCylXZ(:,2) = yCylXZ;


patch('Faces',fCirc,'Vertices',vCylXZ,'FaceColor',[0.85,0.85,0.85],'EdgeColor','none'...);
    ,'specularStrength',specStr,'diffuseStrength',diffStr,'ambientStrength',ambStr,'FaceAlpha',0.7);



%XY plane (red)

xCylXY = sign(yCircV)*500-double(yCircV<0)*500+yCircV;
zCylXY = min(zlim);
yCylXY = zCircV;
vCylXY = [xCylXY,yCylXY, xCylXY];
vCylXY(:,3) = zCylXY;


patch('Faces',fCirc,'Vertices',vCylXY,'FaceColor',[0.85,0.85,0.85],'EdgeColor','none'...);
    ,'specularStrength',specStr/1.8,'diffuseStrength',diffStr/1.8,'ambientStrength',ambStr/1.8,'FaceAlpha',0.7);
%%


hold on;
tStart = timeLagMax+1;
trackColor = [0.3,0.4,0.5];
tStar = tStart;

 
 

    [xS,yS,zS]=sphere(20);
    [Fs,Vs,~]=surf2patch(xS,yS,zS,'Triangles');

    plotHands(1) = patch('Faces',Fs,'Vertices',Vs,'FaceColor',[0.2,0.2,0.2],'EdgeColor','none'...);
    ,'specularStrength',0.2,'diffuseStrength',0.5,'ambientStrength',0.4,'FaceAlpha',0.7);


    alphaData = linspace(0,1,timeLagMax+2); % choose a range of opacities
    %alphaData = linspace(0,1,timeLagMax+2); % choose a range of opacities
    alphaData = alphaData'; % reshape into what matlab expects
    alphaData = (alphaData).^2; % cosmetic options, how 'fast' to fade out

 

   
   plotHands(2) = patchline([xNanPadded(tStar-timeLagMax:tStar)';nan],...
       [yNanPadded(tStar-timeLagMax:tStar)';nan],...
       [zNanPadded(tStar-timeLagMax:tStar)';nan],...
       'edgecolor','k','LineWidth',2.5,...
        'EdgeAlpha','interp',...
        'FaceVertexAlphaData',alphaData);
    
       plotHands(3) = patchline([xNanPadded(tStar-timeLagMax:tStar)';nan],...
       [yNanPadded(tStar-timeLagMax:tStar)';nan],...
       [zNanPadded(tStar-timeLagMax:tStar)';nan],...
       'edgecolor',colorZ,'LineWidth',2.5,...
        'EdgeAlpha','interp',...
        'FaceVertexAlphaData',alphaData);
    
      plotHands(4) = patchline([xNanPadded(tStar-timeLagMax:tStar)';nan],...
       [yNanPadded(tStar-timeLagMax:tStar)';nan],...
       [zNanPadded(tStar-timeLagMax:tStar)';nan],...
       'edgecolor',colorY,'LineWidth',2.5,...
        'EdgeAlpha','interp',...
        'FaceVertexAlphaData',alphaData);

 
      plotHands(5) = patchline([xNanPadded(tStar-timeLagMax:tStar)';nan],...
       [yNanPadded(tStar-timeLagMax:tStar)';nan],...
       [zNanPadded(tStar-timeLagMax:tStar)';nan],...
       'edgecolor',colorX,'LineWidth',2.5,...
        'EdgeAlpha','flat',...
        'FaceVertexAlphaData',alphaData,...
        'FaceAlpha','interp',...
        'FaceColor','none');

 

 
    % loop over remaining timepoints, updating plot object is much faster
    % than recreating
    
    radSp = 60;
    for tStar = timeLagMax+1:size(x,2)+timeLagMax
        newVs = Vs*radSp + [xNanPadded(tStar), yNanPadded(tStar),zNanPadded(tStar)];
        set(plotHands(1),...
            'Vertices',newVs );
          % patchline, which is a line plotter that supports transparency,
          % likes to connect the first and last point together. So
          % including a nan allows one to break this problem and plot just
          % what they want.
        set(plotHands(2),...
            'XData',[xNanPadded(tStar-timeLagMax:tStar)';nan],...
             'YData',[yNanPadded(tStar-timeLagMax:tStar)';nan],...
              'ZData',[zNanPadded(tStar-timeLagMax:tStar)';nan]);
          
        set(plotHands(3),...
            'XData',[xNanPadded(tStar-timeLagMax:tStar)';nan],...
             'YData',[yNanPadded(tStar-timeLagMax:tStar)';nan],...
              'ZData',[(min(zlim)+15)*ones(timeLagMax+1,1);nan]);  
          
          
       set(plotHands(4),...
            'XData',[xNanPadded(tStar-timeLagMax:tStar)';nan],...
             'YData',[(max(ylim)-15)*ones(timeLagMax+1,1);nan],...
              'ZData',[zNanPadded(tStar-timeLagMax:tStar)';nan]);  
          
       set(plotHands(5),...
            'XData',[(min(xlim)+15)*ones(timeLagMax+1,1);nan],...
             'YData',[yNanPadded(tStar-timeLagMax:tStar)';nan],...
              'ZData',[zNanPadded(tStar-timeLagMax:tStar)';nan]);  
          
         figure(gcf);
         %export_fig(gcf,strcat('frame_',num2str(tStar),'.png'))
         pause(0.05); % if grabbing the image of the scene instead of just displaying, this delay can be removed
    end

    

	
	
	
