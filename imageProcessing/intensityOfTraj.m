%traj =csvread('20191213_labeledVsUnlabeledOvernights/glucoseLabeled/20191213_224821.88\5\5_traj_cell_6.csv')
%traj = csvread('ordered50.csv');
%traj =csvread ('20201103_cccp_Everywhere\Cropped\pfvInduced\20201103_234000.24\5\5_traj_cell_9.csv');
%trajFiles = uipickfiles();            
%trajFiles = files2Analyze2(7:end);
%%
%save('C:\Users\dsmendez\Documents\GitHub\shae-miscTools\Diana\20201103_cccp_Everywhere\trajFiles.mat','trajFiles')
%%
%%
%goodTracks= goodTracks(goodIdx);
%dataVars= dataVars(goodIdx);
trajFiles = goodTracks%(goodIdx);
clear meanInt
%trajFiles = posA;
figure
%trajFiles =  trajTresh1000 
%trajFiles =  trajFiles(18)
badIdx=[]
fileIdx = [1:length(goodTracks)]%(goodIdx))]
for iiFiles =fileIdx
    iiFiles
traj = csvread(char(trajFiles(iiFiles)));
idxs=1:length(traj) ;
idxNan =idxs(isnan(traj(:,3)));
clear nextF
for ii=idxNan
    ii
    if ii ==1
        ii2 = 2
        while idxNan(ii2)==ii2
            ii2 = ii2+1
            if ii2>traj(end,6)
                break
            end
        end
        if ii2>traj(end,6)
            break
        end
        prevF=traj(idxs(ii2),3);   
       
    else
        prevF=  traj(idxs(ii)-1,3);
    end
    if ii~=length(traj)
    nextF=  traj(idxs(ii)+1,3);
    end

    iiP=2;
    iiN=2;
    
    while isnan(prevF)
        prevF=  traj(idxs(ii)-iiP,3)
        iiP=iiP+1
    end
    if isempty(nextF)
        nextF =prevF;
    end
    
    while isnan(nextF)
        nextF=  traj(idxs(ii)+iiN,3);
        iiN=iiN+1;
    end
    traj(idxs(ii),3)= (prevF+nextF)/2;
end

try
traj = traj(1:end,:);
x = traj(:,3)*100-(traj(1,3))*100;
%********For newer files******
topInt = traj(:,4);
botInt = traj(:,5);



xo = f.b1;
sigma = f.c1;

A0top = topInt .* exp(((x+deltaZ/2-xo)/sigma).^2);
A0bot = botInt .* exp(((x-deltaZ/2-xo)/sigma).^2);

%%

intArray =[x,A0top, A0bot];

[~,idx] = sort(intArray(:,1)); % sort just the first column
sortedA = intArray(idx,:);   % sort the whole matrix using the sort indices



xCat = horzcat(sortedA(:,1)' ,sortedA(:,1)')';
yCat= horzcat(sortedA(:,2)', sortedA(:,3)')';

%x = horzcat(sortedA(:,1));
%y= horzcat(sortedA(:,3));
edges = (-3000:50:3000);
[~,~,loc]=histcounts(xCat,edges);
meany = accumarray(loc(:),yCat(:))./accumarray(loc(:),1);
xmid = 0.5*(edges(1:end-1)+edges(2:end));
%figure
%plot(x,y,'b.')
hold on; 
plot(xmid(1:size(meany)),meany,'.-')
minInt(iiFiles)=mean(meany(abs(xmid(1:size(meany)))==min(abs(xmid(1:size(meany))))))
meanInt(iiFiles)=mean(yCat);

% 
% figure 
% plot(tempTraj(:,5),x,'p') 

% figure 
% plot(tempTraj(:,5), A0top, '.b')
% hold on
% plot (tempTraj(:,5), A0bot, '.r')
% topInt = traj(:,4);
% botInt = traj(:,5);
% 
% 
% xo = f.b1;
% sigma = f.c1;
% 
% A0top = topInt .* exp(((x+deltaZ/2-xo)/sigma).^2);
% A0bot = botInt .* exp(((x-deltaZ/2-xo)/sigma).^2);
% 
% 
% hold on 
% figure
% plot(x, A0top,\j
% hold on
% plot(x,A0bot,'.r')
catch 
    badIdx = [badIdx,iiFiles]
end
end

%%
figure
histogram(meanInt,5)
 %trajTreshV2000 = trajFiles(meanInt(:)<1500)'
 trajTresh1000 = trajFiles(meanInt(:)<1200)
 goodIdx = meanInt(:)<1200

%%
 trajTresh900 = trajFiles(meanInt(:)<900&meanInt(:)>600)'
 trajTreshfull = trajFiles(meanInt(:)<900)'

%%
figure
plot (x, topInt,'.b')
hold on
plot(x,botInt,'.r')


%%

%traj =csvread('20191213_labeledVsUnlabeledOvernights/glucoseLabeled/20191213_224821.88\5\5_traj_cell_6.csv')
traj = csvread('orderedTest.csv');
x = traj(:,3);
%topInt = traj(:,7);
%botInt = traj(:,8);
partNum = traj(:,5);

for iiPart = 1%:20%:partNum(end)=
    tempTrajFull = traj(partNum ==iiPart,:);
    stIndx=find (tempTrajFull(:,8)==1);
    endIndx =find (tempTrajFull(:,8)==50);
    for counter =1:length(endIndx)
        tempTraj = tempTrajFull(stIndx(counter):endIndx(counter),:);
    
topInt = tempTraj(:,6);
botInt = tempTraj(:,7);

x = tempTraj(:,3)/0.6;
meanx(counter)=mean(x)
xo = f.b1;
sigma = f.c1;

A0top = topInt .* exp(((x+deltaZ/2-xo)/sigma).^2);
A0bot = botInt .* exp(((x-deltaZ/2-xo)/sigma).^2);


hold on 
figure
plot(x, A0top,'.b')
hold on
plot(x,A0bot,'.r')

figure 
plot(tempTraj(:,4),x/0.6,'p') 


    end
end