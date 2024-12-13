%This code corresponds to the calculations and plotting 
%% Load variables 

x20=load('20nmAqLSFixedMSD.mat');
x40=load('40nmPFVFixedMSD.mat');
x50=load('50nmVuldiFixedMSD.mat');
x18= load('40nm-18FixedMSD.mat');
x15= load('40nm+15FixedMSD.mat');
x50exp= load('50nmVuldiExpMSD.mat');
%% Plotting color

color15=1/255*[178,24,43]
color40=1/255*[146,197,222]
color18=1/255*[5,48,97]
color20=1/255*[191,129,45]
color50=1/255*[118,42,131]
color50exp=1/255*[67,24,75]


%% Plots uncorrected MSD
figure 
hold on
plotMSDlog(x15,275,color15,'none')
plotMSDlog(x40,275,color40,'none')
plotMSDlog(x18,275,color18,'none')
plotMSDlog(x20,82,color20,'none')
plotMSDlog(x50,275,color50,'none')
%plotMSDlog(x50exp,color15,'none')
box on
axis square


%% Plot MSD in linear scale

figure 
hold on
plotMSDlinear(x15,275,color15,'none')
plotMSDlinear(x40,275,color40,'none')
plotMSDlinear(x18,275,color18,'none')
%plotMSDlinear(x50exp,275,color50exp,'none')
plotMSDlinear(x20,82,color20,'none')
plotMSDlinear(x50,275,color50,'none')
box on
axis square

%% Fits the MSD data to a diffusion model with dynamic and static localization
%uncertainty
figure
sizeT=60
[a20,b20,c20,rsquare20,sizeV20]=fitLinMSD(x20,15)
%[a20,b20,c20,rsquare20,sizeV20]=fitLinMSD(x20,sizeT)
[a40,b40,c40,rsquare40,sizeV]=fitLinMSD(x40,sizeT)
[a50,b50,c50,rsquare50,sizeV]=fitLinMSD(x50,sizeT)
[a18,b18,c18,rsquare18,sizeV]=fitLinMSD(x18,sizeT)
[a15,b15,c15,rsquare15,sizeV]=fitLinMSD(x15,sizeT)
%[a50exp,b50exp,c50exp,rsquare50exp,sizeV]=fitLinMSD(x50exp,sizeT)


%% Fits MSD with alpha = 0.7
[a20,b20,c20,rsquare20,sizeV20]=fitLinMSD0dot7(x20,sizeT)
[a40,b40,c40,rsquare40,sizeV]=fitLinMSD0dot7(x40,sizeT)
[a50,b50,c50,rsquare50,sizeV]=fitLinMSD0dot7(x50,sizeT)
[a18,b18,c18,rsquare18,sizeV]=fitLinMSD0dot7(x18,sizeT)
[a15,b15,c15,rsquare15,sizeV]=fitLinMSD0dot7(x15,sizeT)


%% Mean diffusion coefficient
figure; 
plot(sizeV20(10:end), a20(10:end)* 10^(-6),'.-g')
hold on
plot(sizeV, a40*10^(-6),'.-r')
plot(sizeV, a50*10^(-6),'.-b') 
plot(sizeV, a18*10^(-6),'.-y')
plot(sizeV, a15*10^(-6),'.-m')
%plot(sizeV, a50exp*10^(-6),'.-k') 
legend('20',...
    '40',...
    '50',...
    '-18',...
    '+15',...
    '50exp')

sigmaSq20=mean(a20(10:13)) 
sigmaSq40=mean(a40(1:end)) 
sigmaSq50=mean(a50(1:end))  
sigmaSq18 =mean(a18(1:end)) 
sigmaSq15 =mean(a15(1:end)) 

%sigmaSq50exp=mean(a50exp(20:40)) 
%%
%% Fitted alpha
figure; plot(sizeV20(10:end), b20(10:end),'.-g')
hold on
plot(sizeV, b40,'.-r')
plot(sizeV, b50,'.-b')
plot(sizeV, b18,'.-y')
plot(sizeV, b15,'.-m')
%plot(sizeV, b50exp,'.-k')
legend('20',...
    '40',...
    '50',...
    '-18',...
    '+15',...
    '50exp')

%% Fitted static localization uncertainty
figure; plot(sizeV20(1:end), c20(1:end),'.-g')
hold on
plot(sizeV, c40,'.-r')
plot(sizeV, c50,'.-b')
plot(sizeV, c18,'.-y')
plot(sizeV, c15,'.-m')
%plot(sizeV, c50exp,'.-k')
legend('20',...
    '40',...
    '50',...
    '-18',...
    '+15',...
    '50exp')

%% Goodness of fit
figure;hold on
%plot(sizeV20(10:end), rsquare20(10:end),'.-g')
%plot(sizeV20(1:end), rsquare20(1:end),'.-g')
plot(sizeV, rsquare40,'.-r')
plot(sizeV, rsquare50,'.-b')
plot(sizeV, rsquare18,'.-y')
plot(sizeV, rsquare15,'.-m')
%plot(sizeV, rsquare50exp,'.-k')
legend('20',...
    '40',...
    '50',...
    '+15',...
    '-18')

%% Average localization uncertainty
sigma40 =mean(c40(5:10))
sigma50=mean(c50(5:10))
sigma18=mean(c18(5:10))
sigma15=mean(c15(5:10))
sigma20=mean(c20(10:13))
%sigma50exp=mean(c50exp(5:10))

sqrt(sigma20/6)
sqrt(sigma40/6)
sqrt(sigma50/6)
%sqrt(sigma50exp/6)
sqrt(sigma15/6)
sqrt(sigma18/6)

%% Fitting 20 nm MSD in x,y and z.
figure
sizeT=60
[a201,b201,c201,rsquare201,sizeV20]=fitLinMSDpartial(x20,20,x20.rowMeanX)
[a202,b202,c202,rsquare202,sizeV20]=fitLinMSDpartial(x20,20,x20.rowMeanY)
[a203,b203,c203,rsquare203,sizeV20]=fitLinMSDpartial(x20,20,x20.rowMeanZ)

sigma201=mean(c201(1:4))
sigma202=mean(c202(1:4))
sigma203=mean(c203(1:4))

endPlot = 7
sqrt(sigma201/2)
sqrt(sigma202/2)
sqrt(sigma203/2)

figure; plot(sizeV20(1:endPlot), c201(1:endPlot),'.-g')
hold on; plot(sizeV20(1:endPlot), c202(1:endPlot),'.-r')
hold on; plot(sizeV20(1:endPlot), c203(1:endPlot),'.-b')

figure; plot(sizeV20(1:endPlot), a201(1:endPlot),'.-g')
hold on; plot(sizeV20(1:endPlot), a202(1:endPlot),'.-r')
hold on; plot(sizeV20(1:endPlot), a203(1:endPlot),'.-b')

figure; plot(sizeV20(1:endPlot), b201(1:endPlot),'.-g')
hold on; plot(sizeV20(1:endPlot), b202(1:endPlot),'.-r')
hold on; plot(sizeV20(1:endPlot), b203(1:endPlot),'.-b')

figure; plot(sizeV20(1:endPlot), rsquare201(1:endPlot),'.-g')
hold on; plot(sizeV20(1:endPlot), rsquare202(1:endPlot),'.-r')
hold on; plot(sizeV20(1:endPlot), rsquare203(1:endPlot),'.-b')
%%
%%Fitting 40 nm -7 MSD in x,y and z.
figure
sizeT=60
[a401,b401,c401,rsquare401,sizeV]=fitLinMSDpartial(x40,sizeT,x40.rowMeanX)
[a402,b402,c402,rsquare402,sizeV]=fitLinMSDpartial(x40,sizeT,x40.rowMeanY)
[a403,b403,c403,rsquare403,sizeV]=fitLinMSDpartial(x40,sizeT,x40.rowMeanZ)

sigma401=mean(c401(1:10))
sigma402=mean(c402(1:10))
sigma403=mean(c403(1:10))

endPlot = 20
sqrt(sigma401/2)
sqrt(sigma402/2)
sqrt(sigma403/2)
figure; plot(sizeV(1:endPlot), c401(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), c402(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), c403(1:endPlot),'.-b')


figure; plot(sizeV(1:endPlot), a401(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), a402(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), a403(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), b401(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), b402(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), b403(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), rsquare401(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), rsquare402(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), rsquare403(1:endPlot),'.-b')


%% Fitting 50 nm MSD in x,y and z.

figure
sizeT=60
[a501,b501,c501,rsquare501,sizeV]=fitLinMSDpartial(x50,sizeT,x50.rowMeanX)
[a502,b502,c502,rsquare502,sizeV]=fitLinMSDpartial(x50,sizeT,x50.rowMeanY)
[a503,b503,c503,rsquare503,sizeV]=fitLinMSDpartial(x50,sizeT,x50.rowMeanZ)

sigma501=mean(c501(3:12))
sigma502=mean(c502(3:12))
sigma503=mean(c503(3:12))

endPlot = 50
sqrt(sigma501/2)
sqrt(sigma502/2)
sqrt(sigma503/2)
figure; plot(sizeV(1:endPlot), c501(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), c502(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), c503(1:endPlot),'.-b')


figure; plot(sizeV(1:endPlot), a501(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), a502(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), a503(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), b501(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), b502(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), b503(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), rsquare501(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), rsquare502(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), rsquare503(1:endPlot),'.-b')

%% Fitting 40 nm +15 MSD in x,y and z.

figure
sizeT=60
[a151,b151,c151,rsquare151,sizeV]=fitLinMSDpartial(x15,sizeT,x15.rowMeanX)
[a152,b152,c152,rsquare152,sizeV]=fitLinMSDpartial(x15,sizeT,x15.rowMeanY)
[a153,b153,c153,rsquare153,sizeV]=fitLinMSDpartial(x15,sizeT,x15.rowMeanZ)

sigma151=mean(c151(3:12))
sigma152=mean(c152(3:12))
sigma153=mean(c153(3:12))

endPlot = 25
sqrt(sigma151/2)
sqrt(sigma152/2)
sqrt(sigma153/2)
figure; plot(sizeV(1:endPlot), c151(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), c152(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), c153(1:endPlot),'.-b')


figure; plot(sizeV(1:endPlot), a151(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), a152(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), a153(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), b151(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), b152(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), b153(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), rsquare151(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), rsquare152(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), rsquare153(1:endPlot),'.-b')


%% Fitting 40 nm -18 MSD in x,y and z.

figure
sizeT=60
[a181,b181,c181,rsquare181,sizeV]=fitLinMSDpartial(x18,sizeT,x18.rowMeanX)
[a182,b182,c182,rsquare182,sizeV]=fitLinMSDpartial(x18,sizeT,x18.rowMeanY)
[a183,b183,c183,rsquare183,sizeV]=fitLinMSDpartial(x18,sizeT,x18.rowMeanZ)

sigma181=mean(c181(1:20))
sigma182=mean(c182(1:20))
sigma183=mean(c183(1:20))

endPlot = 20
sqrt(sigma181/2)
sqrt(sigma182/2)
sqrt(sigma183/2)
figure; plot(sizeV(1:endPlot), c181(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), c182(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), c183(1:endPlot),'.-b')


figure; plot(sizeV(1:endPlot), a181(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), a182(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), a183(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), b181(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), b182(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), b183(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), rsquare181(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), rsquare182(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), rsquare183(1:endPlot),'.-b')


%% Fitting 50nm exp MSD in x,y and z.

figure
sizeT=60
[a50exp1,b50exp1,c50exp1,rsquare50exp1,sizeV]=fitLinMSDpartial(x50exp,sizeT,x50exp.rowMeanX)
[a50exp2,b50exp2,c50exp2,rsquare50exp2,sizeV]=fitLinMSDpartial(x50exp,sizeT,x50exp.rowMeanY)
[a50exp3,b50exp3,c50exp3,rsquare50exp3,sizeV]=fitLinMSDpartial(x50exp,sizeT,x50exp.rowMeanZ)


[a50exp3,b50exp3,c50exp3,rsquare50exp3,sizeV]=fitLinMSD0dot7partial(x50exp,sizeT,x50exp.rowMeanZ)


sigma50exp1=mean(c50exp1(1:20))
sigma50exp2=mean(c50exp2(1:20))
sigma50exp3=mean(c50exp3(1:20))

endPlot = 30
sqrt(sigma50exp1/2)
sqrt(sigma50exp2/2)
sqrt(sigma50exp3/2)
figure; plot(sizeV(1:endPlot), c50exp1(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), c50exp2(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), c50exp3(1:endPlot),'.-b')


figure; plot(sizeV(1:endPlot), a50exp1(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), a50exp2(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), a50exp3(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), b50exp1(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), b50exp2(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), b50exp3(1:endPlot),'.-b')

figure; plot(sizeV(1:endPlot), rsquare50exp1(1:endPlot),'.-g')
hold on; plot(sizeV(1:endPlot), rsquare50exp2(1:endPlot),'.-r')
hold on; plot(sizeV(1:endPlot), rsquare50exp3(1:endPlot),'.-b')




%% xy fitting 
figure
sizeT=60
%[a50expxy,b50expxy,c50expxy,rsquare50expxy,sizeV]=fitLinMSD0dot7partial(x50exp,sizeT,x50exp.rowMeanX+x50exp.rowMeanY-sigmaxy50exp)
[a50xy,b50xy,c50xy,rsquare50xy,sizeV]=fitLinMSD0dot7partial(x50,sizeT,x50.rowMeanX+x50.rowMeanY-sigmaxy50)
[a40xy,b40xy,c40xy,rsquare40xy,sizeV]=fitLinMSD0dot7partial(x40,sizeT,x40.rowMeanX+x40.rowMeanY-sigmaxy40)
[a20xy,b20xy,c20xy,rsquare20xy,sizeV]=fitLinMSD0dot7partial(x20,sizeT,x20.rowMeanX+x20.rowMeanY-sigmaxy20)
[a18xy,b18xy,c18xy,rsquare18xy,sizeV]=fitLinMSD0dot7partial(x18,sizeT,x18.rowMeanX+x18.rowMeanY-sigmaxy18)
[a15xy,b15xy,c15xy,rsquare15xy,sizeV]=fitLinMSD0dot7partial(x15,sizeT,x15.rowMeanX+x15.rowMeanY-sigmaxy15)



%% 
%% Plotting 2D MSDs 
sigmaxy20 = sigma201+sigma202;
sigmaxy40 = sigma401+sigma402;
sigmaxy50 = sigma501+sigma502;
%sigmaxy50exp = sigma50exp1+sigma50exp2;
sigmaxy15 = sigma151+sigma152;
sigmaxy18 = sigma181+sigma182;
figure 
hold on
box on

%D50exp=mean(a50expxy(7:13))*10^-6
D50=mean(a50xy)*10^-6
D40=mean(a40xy(2:13))*10^-6
D18=mean(a18xy)*10^-6
D15=mean(a15xy)*10^-6


color15=1/255*[178,24,43]
plotMSDlinear2Dcorr(x15,275,color15,sigmaxy15-0.0035*10^6,'-')
color40=1/255*[146,197,222]
plotMSDlinear2Dcorr(x40,275,color40,sigmaxy40-0.0035*10^6,'-')
color18=1/255*[5,48,97]
plotMSDlinear2Dcorr(x18,275,color18,sigmaxy18-0.0035*10^6,'-')
color20=1/255*[191,129,45]
 plotMSDlinear2Dcorr(x20,125,color20,sigmaxy20-0.0035*10^6,'-')
 color50=1/255*[118,42,131]
 plotMSDlinear2Dcorr(x50,350,color50,sigmaxy50-0.0035*10^6,'-')
% color50exp=1/255*[67,24,75]

set(gca,'xscale','lin')
set(gca,'yscale','lin')


%% Corrected MSD log plot
figure
plotMSDlogCorr(sigma20,x20,82,color20,1)
hold on
plotMSDlogCorr(sigma40,x40,275,color40,1)
plotMSDlogCorr(sigma50,x50,275,color50,1)
%plotMSDlogCorr(sigma50exp,x50exp,275,color50exp,1)
 plotMSDlogCorr(sigma18,x18,275,color18,1)
 plotMSDlogCorr(sigma15,x15,275,color15,1)
box on
axis square
xlabel('\tau (s)')
ylabel('MSD (\mum^2)')


%% Corrected localization error and normalized by diffusion coeff
figure
sizeT=275
idxPlot =5
color=1/255*[146,197,222]
kCycle = x40.kCycle;
rowMeanFull=(x40.rowMeanFull(idxPlot:end)-sigma40);
xlin=(idxPlot:idxPlot+sizeT-1)*kCycle

figure
plot(xlin, rowMeanFull(1:sizeT)/(sigmaSq40), 'LineStyle','-',...
    'LineWidth',3,'Color',color)
set(gca, 'FontSize', 20)%,'BinWidth', 50)


color=1/255*[118,42,131]
kCycle = x50.kCycle;
rowMeanFull=(x50.rowMeanFull(idxPlot:end)-sigma50);

hold on
plot(xlin, rowMeanFull(1:sizeT)/(sigmaSq50), 'LineStyle','-',...
    'LineWidth',3,'Color',color)
set(gca, 'FontSize', 20)%,'BinWidth', 50)


%color='k'
%kCycle = x50exp.kCycle;
%rowMeanFull=(x50exp.rowMeanFull(idxPlot:end)-sigma50exp);

%hold on
%plot(xlin, rowMeanFull(1:sizeT)/(sigmaSq50exp), 'LineStyle','-',...
%    'LineWidth',3,'Color',color)
%set(gca, 'FontSize', 20)%,'BinWidth', 50)



color=1/255*[5,48,97]
kCycle = x18.kCycle;

rowMeanFull=(x18.rowMeanFull(idxPlot:end)-sigma18);

hold on
plot(xlin, rowMeanFull(1:sizeT)/(sigmaSq18), 'LineStyle','-',...
    'LineWidth',3,'Color',color)
set(gca, 'FontSize', 20)%,'BinWidth', 50)


color=1/255*[178,24,43]
kCycle = x15.kCycle;
rowMeanFull=(x15.rowMeanFull(idxPlot:end)-sigma15);



hold on
plot(xlin, rowMeanFull(1:sizeT)/(sigmaSq15), 'LineStyle','-',...
    'LineWidth',3,'Color',color)
set(gca, 'FontSize', 20)%,'BinWidth', 50)

sizeT =87
idxPlot= 5
color=1/255*[191,129,45]
kCycle = x20.kCycle;
rowMeanFull=(x20.rowMeanFull(idxPlot:end)-sigma20);
xlin=(idxPlot:idxPlot+sizeT-1)*kCycle

hold on
plot(xlin, rowMeanFull(1:sizeT)/(sigmaSq20), 'LineStyle','-',...
    'LineWidth',3,'Color',color)
set(gca, 'FontSize', 20)%,'BinWidth', 50)
 set(gca,'xscale','log')
 set(gca,'yscale','log')



%% Removing the bad trajectories of x18
  x18.goodTracks(2)=[]
  x18.goodTracks(6)=[]
  x18.goodTracks(6)=[]
  
  
  x18.dataVars(2)=[]
  x18.dataVars(6)=[]
  x18.dataVars(6)=[]
%% Find the particles in the poles and nucleoid
[x18Out,~,x18Pole]=findOutNucNoNan(x18);
[x40Out,~,x40Pole]=findOutNucNoNan(x40);
[x50Out,~,x50Pole]=findOutNucNoNan(x50);
[x15Out,~,x15Pole]=findOutNucNoNan(x15);
[x20Out,~,x20Pole]=findOutNucNoNan(x20);

%% Calculate the fraction of particles in the pole
x20PFrac=(nansum(x20Pole(:))/length(x20Pole(:)))
x40PFrac=(nansum(x40Pole(:))/length(x40Pole(:)))
x50PFrac=(nansum(x50Pole(:))/length(x50Pole(:)))
x18PFrac=(nansum(x18Pole(:))/length(x18Pole(:)))
x15PFrac=(nansum(x15Pole(:))/length(x15Pole(:)))

%% Calculate the fraction of particles in nucleoid
x20Frac=1-(nansum(x20Out(:))/length(x20Out(:)))
x40Frac=1-(nansum(x40Out(:))/length(x40Out(:)))
x50Frac=1-(nansum(x50Out(:))/length(x50Out(:)))
x18Frac=1-(nansum(x18Out(:))/length(x18Out(:)))
x15Frac=1-(nansum(x15Out(:))/length(x15Out(:)))


%% Calculate the standard deviation of fraction in nucleoid using bootstrap
intIn20 =bootci(200,{@(x)1-(nansum(x)/length(x)),x20Out(:)'})';
intIn40 =bootci(200,{@(x)1-(nansum(x)/length(x)),x40Out(:)'})';
intIn50 =bootci(200,{@(x)1-(nansum(x)/length(x)),x50Out(:)'})';
intIn180 =bootci(200,{@(x)1-(nansum(x)/length(x)),x18Out(:)'})';
intIn150 =bootci(200,{@(x)1-(nansum(x)/length(x)),x15Out(:)'})';


%% Calculate the standard deviation of fraction in poles using bootstrap
inPoleSD20 =bootci(200,{@(x)(nansum(x)/length(x)),x20Pole(:)'})';
inPoleSD40 =bootci(200,{@(x)(nansum(x)/length(x)),x40Pole(:)'})';
inPoleSD50 =bootci(200,{@(x)(nansum(x)/length(x)),x50Pole(:)'})';
inPoleSD180 =bootci(200,{@(x)(nansum(x)/length(x)),x18Pole(:)'})';
inPoleSD150 =bootci(200,{@(x)(nansum(x)/length(x)),x15Pole(:)'})';


%% Plot pole occupation figure
X1 = [20,40,50]
Y2 = [x20PFrac*100,x40PFrac*100,x50PFrac*100]
y2min=[inPoleSD20(1),inPoleSD40(1),inPoleSD50(1)]*100
y2max=[inPoleSD20(2),inPoleSD40(2),inPoleSD50(2)]*100
lw=0.8

createInPoleFig(X1, Y2, y2min, y2max,lw,{'20','30','40','50'},'Diameter (nm)',[20,30,40,50],color20,color40,color50)

X1 = [-18,-7,+15]
Y2 = [x18PFrac*100,x40PFrac*100,x15PFrac*100]
y2min=[inPoleSD180(1),inPoleSD40(1),inPoleSD150(1)]*100
y2max=[inPoleSD180(2),inPoleSD40(2),inPoleSD150(2)]*100


createInPoleFig(X1, Y2, y2min, y2max,lw,{'-2160','-840','0','+1800'},'Charge (e)',[-18,-7,0,+15],color18,color40,color15)


