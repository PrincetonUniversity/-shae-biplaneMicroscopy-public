function [a,b,c,rsquare,sizeTvec]=fitLinMSD0dot7(struct2Analize,sizeT)
index = 1;
clear rsquare a b c 
sizeTvec =3:sizeT;
for ii=sizeTvec
    try
    x = struct2Analize.kCycle*(1:ii)';
    y = struct2Analize.rowMeanFull(1:ii);%*10^-6;
    [temp1,temp2]=  fitMSD0dot7(x, y)
    a(index)=temp1.a; 
    b(index)=1;
    c(index)=temp1.c ;
    gofFit(index)= temp2;
    rsquare(index)=temp2.rsquare;
    index=index+1;
    catch
        a(index)=NaN;
        b(index)=NaN;
        c(index)=NaN;
        gofFit(index)=gofFit(index-1);
        rsquare(index)=NaN;
    end
end