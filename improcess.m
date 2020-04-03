%AUTHOR:
    %DANIEL TOVBIS (2019)
%DESCRIPTION:
    %This function does preprocessing for the image defined by imname. It is
    %mostly called by other functions in the pipeline and does not typically need to be run
    %separately. Using regist4x75 should have all these values already populated.
    %Steps include an increase in local contrast, sharpening, and conversion to
    %grayscale. Then diffusion filtering, and opening-and-closing by
    %reconstruction.
%INPUTS:
    %imname: A workspace variable containing the image to process
    %erodesize: Size of structuring element used for morphological erosion.
    %dilatesize: Size of structuring elment used for morphological dilation.
%OUTPUTS:
    %Iobrcbr: The processed image.
function Iobrcbr=improcess(imname,erodesize,dilatesize)
currentdir=pwd;
rgb= imname;
if size(rgb,3)==3 %if the image is rgb (which it normally shouldn be)
    rgblc=localcontrast(rgb,0.3,0.6);
    rgbs=imsharpen(rgblc);
    I = rgb2gray(rgbs);
else
    I=imname; %if the image is already grayscale (which it normally shouldn't be)
end
%imshow(I)
%% Image preprocessing
I=imdiffusefilt(I,'NumberofIterations',5,'Conductionmethod','exponential'); %diffusion filter
%Create structuring element
se = strel('disk',erodesize,0); %set up structural element
%Image eroding
Ie = imerode(I, se); %erode
%Reconstruct image using erosion to enhance contrast
Iobr = imreconstruct(Ie, I);
se2=strel('disk',dilatesize,0); %set up structural element
Iobrd = imdilate(Iobr, se2); %dilate
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr)); %reconstruct using image complements
Iobrcbr = imcomplement(Iobrcbr); %invert
end