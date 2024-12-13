%attempt at fitting a cell 

function [F,frot] = fitCellonlyEndCapsNoConvRot(a,data)

%Inputs:
% a = vector of parameters. 
%     a(1)= l
%     a(2)= r
%     a(3)= x0 center cell 
%     a(4)= y0 center cell 
%     a(5)= normIntensity
%     a(6)= theta(for end caps) 
%
% data = 2d matrix of image data or 2d matrix which is the same size as
%        image data
       

%% functional form of 2d gaussian from Wikipedia
if (abs(a(7))>45)
    data=data';
end
sizex = size(data,2);
sizey = size(data,1);

%% cell shape
l=a(1);
r=a(2);
x0 =a(3);
y0=a(4);
th =a(6);
%% Drawing the midcell
if (abs(a(7))>45)
    y0 =a(3)
    x0=a(4)
end

frec =drawRectangle(sizex,sizey,l,2*r,x0,y0);

% [x,y]= meshgrid(1:sizex,1:sizey);
% fx= (x-x0>-l/2 & x-x0 <l/2);
% 
% fy=(y-y0>-w/2& y-y0<w/2);
% rect= fx.*fy;




%% Drawing the caps
cr=x0+l/2;
cl=x0-l/2;

d = r/tand(th);
Rc = r/sind(th);

xcr = cr-d ;
xcl=cl+d;
yc=y0;

circMaskr = drawCircle(sizex,sizey,Rc,xcr,yc);
circMaskl = drawCircle(sizex,sizey,Rc,xcl,yc);


recMask = drawRectangle(sizex,sizey,l,sizey,x0,y0);

fcapR = circMaskr .* ~recMask;
fcapL = circMaskl .* ~recMask;

ftotal = fcapR +frec+ fcapL;

%%
%ftotal2 =makeCellMask(sizex,sizey,r, l);

%ffilt =imgaussfilt(ftotal,2);

if (abs(a(7))>45)
    ftotal=ftotal';
%     bwmem =bwperim(ftotal);
%     se90 = strel('line',1,90);
%     se0 = strel('line',1,0);
%     bwInMemDil = imdilate(bwmem,[se90 se0]); %close edges
%     figure; imshow(bwInMemDil,[])
    if a(7)>0
        frot = imrotate(ftotal,a(7)-90,'bilinear','crop');
    else 
        frot = imrotate(ftotal,a(7)+90,'bilinear','crop');
    end
else
%    bwmem =bwperim(ftotal);
%     se90 = strel('line',1,90);
%     se0 = strel('line',1,0);
%     bwInMemDil = imdilate(bwmem,[se90 se0]); %close edges
%     figure; imshow(bwInMemDil,[])
    frot = imrotate(ftotal,a(7),'bilinear','crop');
end


% frot = imrotate(ftotal,a(7),'bilinear','crop');

fcell =imgaussfilt(frot,2);
%a(5)
5
F =a(5)*fcell;
imshow(F,[])
drawnow


