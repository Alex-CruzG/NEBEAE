function [ classificationColorMap ] = generateColorMap( classificationMap )
%GENERATECOLORMAP 
%   This code generate a color map based in the results obtained from the
%   HELICoiD classification result. It assigns colors depending on a class
%   list and the label in each pixel of the classification map.
%   
%   classificationMap:  2D Matrix with the labeled classification results.
%
%
% Other m-files required: none
% Subfunctions: classList.m
% MAT-files required:
%
% Authors: Samuel Ortega, Himar Fabelo
% email address: sortega@iuma.ulpgc.es, hfabelo@iuma.ulpgc.es
% November 2017.

    sizeHC = size(classificationMap);
    classificationColorMap = zeros(sizeHC(1), sizeHC(2), 3, 'uint8');
    classList = getLabelList;
    
    for i = 1 : sizeHC(1)
        for j = 1 : sizeHC(2)
            currentClassIndex = strcmp(num2str( classificationMap( i, j ) ), classList( :, 1 ) );
            currentClassColor = classList( currentClassIndex, 5 );
            classificationColorMap( i, j, : ) = currentClassColor{ 1 };
        end
    end

end

%% 

function [ sampleClassName ] = getLabelList()
%SAMPLECLASSES Summary of this function goes here
% This function generate the list of possible materials/substances
% that can be presented in the image
    
    sampleClassName = {
    '0','Noinfo','Noinfo','Noinfo',[255 255 255];
    '1','Normal','NonDefined','NonDefined',[0 255 0];
    '2','Tumour','Primary-GIV','PureGBM',[255 0 0];
    '3','Other','Blood','Generic',[0 0 255]; 
    '4','Background','BG-Generic','BG-Generic',[0 0 0];
 };

end