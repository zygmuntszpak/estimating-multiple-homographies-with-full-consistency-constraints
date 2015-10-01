%   Function: compute_mean_root_mean_square_error
%
%   The root-mean-square error is computed for residuals associated 
%   with each homography  in listOfErrors and then finally, a further
%   average across all homographies is computed. 
%
%   Parameters:
%
%      listOfErrors           - a data structure containing squared
%                               residuals. The first index designates
%                               the scene, and the second index designates
%                               a particular homography
%
%
%
%   Returns:
%
%    A scalar representing the mean root-mean-square error computed 
%    by averaging the errors over all scenes and all planes (homographies), 
%    and a matrix containing the root-means-square errors for each
%    plane by averaging over the number of scenes 
%
%   See Also:
%
%  visualise_symmetric_transfer_errors_for_fixed_noise
%
%
%  Zygmunt L. Szpak (c) 2014
%  Last modified 16/09/2014
function [ meanRmsError, rmsErrorForEachPlane] = ...
                        compute_mean_root_mean_square_error(listOfErrors)

    [numberOfScenes, numberOfPlanes] = size(listOfErrors);
    
    meanSquaredErrorForEachPlane = zeros(numberOfScenes,numberOfPlanes);
        
    for i=1:numberOfScenes
        for j=1:numberOfPlanes 
        % Since each point is represented by an x and y coordinate we
        % divide the total number of entries in the vector by 2 
        numberOfPoints = length(listOfErrors{i,j}.errorsForEachPlane) / 2;
        % we divide by four since we average over the (x,y) coordinates
        % in the first image, and the (x,y) coordinates in the second image
        meanSquaredErrorForEachPlane(i,j) = ...
              sum(listOfErrors{i,j}.errorsForEachPlane)/(numberOfPoints*4);        
        end
    end
    
    % average over the number of scenes
    meanSquaredErrorForEachPlane = mean(meanSquaredErrorForEachPlane);
    
    % take the square root for each plane
    rmsErrorForEachPlane = sqrt(meanSquaredErrorForEachPlane);
    
    % now take the average over all the planes
    meanRmsError = mean(rmsErrorForEachPlane);
    
    

end

