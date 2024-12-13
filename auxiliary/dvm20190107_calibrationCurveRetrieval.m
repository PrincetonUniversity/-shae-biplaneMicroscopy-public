%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [nameRed,pathRed]=uigetfile('*.csv','Choose red file');
% redCSV=fullfile(pathRed,nameRed);
% [nameBlue,pathBlue]=uigetfile('*.csv','Choose blue file');
% blueCSV=fullfile(pathBlue,nameBlue);
%%
FRAMES =40% 36;%360;
%%
redFiles =uipickfiles';
blueFiles =uipickfiles';
%%
for iiFiles =1:size(redFiles,1)
    if size(redFiles,1) ~=size(blueFiles,1)
        break;
    end

    meanBGRed= csvread(char(redFiles(iiFiles)));
    meanBGBlue= csvread(char(blueFiles(iiFiles)));
   
     nanIndx = find(isnan(meanBGRed)==1);
    [rowNan,colNan]=ind2sub((size(meanBGRed)),nanIndx);
    meanBGRed(rowNan,:)=[];
    meanBGBlue(rowNan,:)=[];
    nanIndx = find(isnan(meanBGBlue)==1);
    [rowNan,colNan]=ind2sub((size(meanBGRed)),nanIndx);
    meanBGRed(rowNan,:)=[];
    meanBGBlue(rowNan,:)=[];
     
    
    normBGRed = meanBGRed./repmat(max(meanBGRed,[],2),1,FRAMES);%30);
    normBGBlue = meanBGBlue./repmat(max(meanBGBlue,[],2),1,FRAMES);%30);

        clear minV1
        clear minV2
        clear shiftedRed
        clear shiftedBlue
        for iiImage = 1:size(meanBGRed,1)
            min1 = 100000000000000;
            min2 = 100000000000000;
            for iiDisp = -40:40
                temp1 =sum (sqrt((normBGRed(1,:) -circshift(normBGRed(iiImage,:),iiDisp)).^2));
                temp2 =sum (sqrt((normBGBlue(1,:) -circshift(normBGBlue(iiImage,:),iiDisp)).^2));
                if temp1 < min1 && temp2<min2
                    min1 =temp1;
                    min2= temp2;
                    disp1 = iiDisp;
                end
            end
            minV1(iiImage)=min1;
            minV2(iiImage)=min2;
            optDips(iiImage) =disp1;
         end



    %%
    outliers =  minV1+minV2>quantile(minV1+minV2,0.95);
    normBGRed(minV1+minV2>quantile(minV1+minV2,0.95),:)=[];
    normBGBlue(minV1+minV2>quantile(minV1+minV2,0.95),:)=[];
    optDips(:,minV1+minV2>quantile(minV1+minV2,0.95))=[];


    %%
    for iiBead = 1:size(normBGRed,1)
        shiftedRed(iiBead,:) = circshift(normBGRed(iiBead,:),optDips(iiBead));
        shiftedBlue(iiBead,:) = circshift(normBGBlue(iiBead,:),optDips(iiBead));
    end


%    XRed = -1790:10:1800;
%    XBlue = -1790:10:1800;

   XRed = -1700:100:1800;
   XBlue = -1700:100:1800;

    meanBlue=(mean(shiftedBlue));
    meanRed=(mean(shiftedRed));
    outIndxBlue=isnan(meanBlue);


    meanBlue(outIndxBlue)=[];
    meanRed(outIndxBlue)=[];
    XRed(outIndxBlue)=[];
    XBlue(outIndxBlue)=[];
    shiftedBlue(:,outIndxBlue)=[];
    shiftedRed(:,outIndxBlue)=[];

    outIndxRed=isnan(meanRed);
    meanBlue(outIndxRed)=[];
    meanRed(outIndxRed)=[];
    XRed(outIndxRed)=[];
    XBlue(outIndxRed)=[];
    shiftedBlue(:,outIndxRed)=[];
    shiftedRed(:,outIndxRed)=[];
figure
plot(meanBlue,'.')
hold on
plot(meanRed,'.')

    %%
    % 
    % limitRangeIndx=(XRed<400&XRed>-550);
    % 
    % RedBeads = shiftedRed(:,limitRangeIndx);
    % BlueBeads= shiftedBlue(:,limitRangeIndx);
    % 
    % logDiv=log(meanBlue./meanRed);
    % figure
    % hold on;
%     plot(logDiv(limitRangeIndx),XBlue(limitRangeIndx),'pb');
%     
    % 
    % 
    % ft = fittype( 'cubicFit( x, a, b, c, d )' );
    % f = fit( logDiv(limitRangeIndx)', XBlue(limitRangeIndx)', ft, 'StartPoint', [-1.5*10^2, -2.1*10^2, 24, -10] );
    % plot( f) 
    % save('calibrationCurve','f')
    % 
    % figure
    % hold on;
    % plot(XBlue(limitRangeIndx),f(logDiv(limitRangeIndx))-XBlue(limitRangeIndx)')
    % mean(f(logDiv(limitRangeIndx))-XBlue(limitRangeIndx)')
    % std(f(logDiv(limitRangeIndx))-XBlue(limitRangeIndx)')
    % %%
    % 
    % limitRangeIndx=(XRed<400&XRed>-550);
    % 
    % RedBeads = shiftedRed(:,limitRangeIndx);
    % BlueBeads= shiftedBlue(:,limitRangeIndx);
    % 
    % logDiv=log(meanBlue./meanRed);
    % figure
    % hold on;
    % plot(logDiv(limitRangeIndx),XBlue(limitRangeIndx),'pb');
    % 
    % 
    % 
    % ft = fittype( 'cuarticFit( x, a, b, c, d,e )' );
    % f = fit( logDiv(limitRangeIndx)', XBlue(limitRangeIndx)', ft, 'StartPoint', [-1.5*10^2, -2.1*10^2, 24, -10,0.37] );
    % plot( f) 
    % save('calibrationCurve','f')
    % 
    % figure
    % plot(XBlue(limitRangeIndx),f(logDiv(limitRangeIndx))-XBlue(limitRangeIndx)')
    % mean(f(logDiv(limitRangeIndx))-XBlue(limitRangeIndx)')
    % std(f(logDiv(limitRangeIndx))-XBlue(limitRangeIndx)')

    %%
    %limitRangeIndx=(XRed<800&XRed>-300); %Voll5
    %  limitRangeIndx=(XRed<701&XRed>-101); %Full Range
   %limitRangeIndx=(XRed<1801&XRed>-1701); %Full Range
    limitRangeIndx=(XRed<401&XRed>-501); %voll the one
    %limitRangeIndx=(XRed<1001&XRed>-201); %good for intensity
    RedBeads = shiftedRed(:,limitRangeIndx);
    BlueBeads= shiftedBlue(:,limitRangeIndx);
    
    
    negIndx = find(RedBeads<=0);
    [rowNeg,colNeg]=ind2sub((size(RedBeads)),negIndx);
    RedBeads(rowNeg,:)=[];
    BlueBeads(rowNeg,:)=[];
    negIndx = find(BlueBeads<=0);
    [rowNeg,colNeg]=ind2sub((size(BlueBeads)),negIndx);
    RedBeads(rowNeg,:)=[];
    BlueBeads(rowNeg,:)=[];
    

    logDivision = log(RedBeads./BlueBeads);
    logDiv=log(meanBlue./meanRed);

     infIndx = find(isinf(logDivision)==1);
    [rowInf,colInf]=ind2sub((size(logDivision)),infIndx);
    logDivision(rowInf,:)=[];
    nanIndx = find(isnan(logDivision)==1);
    [rowNan,colNan]=ind2sub((size(logDivision)),infIndx);
    logDivision(rowNan,:)=[];

    figure
    hold on;
    plot(mean(logDivision),XBlue(limitRangeIndx),'pb');



    ft = fittype( 'grade5Fit( x, a, b, c, d,e,f )' );
    f = fit( mean(logDivision)', XBlue(limitRangeIndx)', ft, 'StartPoint', [-1.5*10^2, -2.1*10^2, 24, -10,0.37,1.1] );
    plot( f) 
    save('calibrationCurveVoll5','f')

    %Testing the fit
   figure
   plot(XBlue(limitRangeIndx),f(mean(logDivision))-XBlue(limitRangeIndx)')
    mean(f(mean(logDivision))-XBlue(limitRangeIndx)');
    std(f(mean(logDivision))-XBlue(limitRangeIndx)');

    %% Brenner gradient ratio plot     

    figure
   hold on;
    logDivision = log(RedBeads./BlueBeads);
    E=std(logDivision);
   errorbar(XBlue(limitRangeIndx),mean(logDivision),E);
   ylabel('log ratio');
   xlabel('10 nm steps (after focal shift)');
    plot(XBlue(limitRangeIndx),mean(logDivision),'pb');
    
    
    
    % %% Brenner gradient ratio plot     
    % 
    % figure
    % hold on;
    % logDivision= log(RedBeads./BlueBeads);
    % for ii=1:size(logDivision,2)
    %     logDivisionPrime(:,ii) =f( logDivision(:,ii));
    % end
    % E=std(logDivisionPrime);
    % errorbar(XBlue(limitRangeIndx),mean(logDivisionPrime),E);
    % ylabel('log ratio');
    % xlabel('10 nm steps (after focal shift)');


    maxR= max(mean(logDivision));
    minR=min(mean(logDivision));
    range= maxR-minR;
    sum(E)/range;
   
%     hold on
%     plot(XBlue(limitRangeIndx),E/range,'-r')
%     totalDev(iiFiles)=sum(E)/range;

end