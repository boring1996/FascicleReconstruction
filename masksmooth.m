%AUTHOR:
    %JACK HAN (2019)
%DESCRIPTION:
    %This function smooths a binary mask using a convultional mean filter 
    %by converting it to a single precision array, then rethresholds to
    %recreate a binary image.
%INPUTS:
    %binaryImage: the mask to smooth
    %windowSize: the size of the filter
    %thresh: The threshold to use after smoothing to recreate a binary
    %image.
%OUTPUTS:
    %smoothedImage: A smoothed binary mask.
function smoothedImage = masksmooth(binaryImage,windowSize,thres)
kernel = ones(windowSize) / windowSize ^ 2; %Create kernel
blurryImage = conv2(single(binaryImage), kernel, 'same'); %Apply mean filter
smoothedImage = blurryImage > thres; % Rethreshold