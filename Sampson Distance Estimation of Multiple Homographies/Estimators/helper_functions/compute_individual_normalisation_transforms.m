function [listOfT listOfTp] = compute_individual_normalisation_transforms(sceneData)
    numOfH = length(sceneData);
    listOfT = cell(1,numOfH);
    listOfTp = cell(1,numOfH);  
    
    for i = 1:numOfH
        x = sceneData(i).ptsInViewOne;
        xp = sceneData(i).ptsInViewTwo;
        n = length(x);
        x = [ x; ones( 1, n ) ];
        xp = [ xp; ones( 1, n ) ];
        T = normalize_data_transform( x );
        Tp = normalize_data_transform( xp );
        listOfT{i} = T;
        listOfTp{i} = Tp;
    end

end