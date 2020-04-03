%AUTHOR:
    %DANIEL TOVBIS (2019)
%DESCRIPTION:
    %Simple Preprocessing for iHC images, including otsu thresholding, hole filling, image closing,
    %and area thresholding.
%INPUTS:
    %imname: Workspace variable containing the image to process.
    %Populated by default as part of regist4x_75.
%OUTPUTS:
    %contour: The processed image.
function contour=improcessihc(imname)
currentdir=pwd;
rgb= imname;
end
if size(rgb,3)==3
    I = rgb2gray(rgb);
else
    I=imname;
end
se=strel('disk',3);
se2=strel('disk',5);
[counts,~]=imhist(I,64);
T=otsuthresh(counts);
bw1=imcomplement(imbinarize(I,T));
bw15=imfill(bw1,'holes'); %Fill holes
bw2=imclose(bw15,se); %Morph closing
bw25=imfill(bw2,'holes'); %File holes again
bw3=imopen(bw25,se2); %Opening
bw4=bwareaopen(bw3,250); %Remove any small components
contour=bw4;
end
