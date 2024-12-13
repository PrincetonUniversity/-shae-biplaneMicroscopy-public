function tiffListFull = chooseTiffsMulticolor(expDir)
   
    subExp = getDirectoryNames(expDir)
    5;
    
    tiffListFull =[];
    for iiExp= 1:numel(subExp)
        zStackPath = strcat (expDir, '/',char(subExp(iiExp)));
        imageType = getDirectoryNames(zStackPath);

        a=strfind(imageType,'GEM');
        for iiStep =1:length(a)
            if ~isempty(a{iiStep})
                fullPathTif= cell2mat(fullfile(expDir, subExp(iiExp),imageType(iiStep),'*.tif'));
                fullPath= cell2mat(fullfile(expDir, subExp(iiExp),imageType(iiStep)));
                dd = dir(fullPathTif);
                tifList = {dd.name};
                tifList = fullfile(fullPath,tifList);
                tiffListFull = [tiffListFull, tifList]
                tiffListFull'
            end
        end

    end

end

