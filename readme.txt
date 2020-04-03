SETUP:
Your files should be organized as follows:

Code/ImageDirectory/Ground Truths

All scripts and functions in Code.
ImageDirectory has folders with images separated by specimen (as in, Distal1IHC, Proximal1HE, etc)
Ground Truths contains the ground truth segmentations.
ImageDirectory also contains training data.

SCRIPTS:
TheGrandImageRegistrationScript contains pre-set parameters for registering all images
CreateMaskforallimages_allmethods contains pre-set parameters for segmentation/detection for all images.
fascrecon_modelonly generates the models.

FUNCTIONS:
improcess/improcessihc for pre-processing H&E and IHC images
regist4x75 for registration
NN_ACoutsidein for detection/segmentation (H&E)
masksmooth smooths masks after resizing.
interpmask does interpolation.
tranFascicleRCNN_Noagument_NSG trains the neural networks using the training data given.

DATA STRUCTURE:
The training data is set up to find the images using a preset directory.
So you will need to use strrep to change the data source of the groundTruth object to whatever directory structure you have the images saved in.
EG:
currentdir=pwd
N=size(gTruthDistal3HE.DataSource,1);
for q=1:N
    gTruthDistal3HE.DataSource{q}=strrep(gTruthDistal3HE.DataSource{q},'M:\Peripheral Nerve Studies\Daniel\Nerve Processing',currentdir);
end

This code will change the M:\Peripheral Nerve Studies\Daniel\Nerve Processing\Nerve Images to currentdir\Nerve Images\, where currentdir is your current working directory
Code to make these changes has been included in the TrainFascicleRCNN script.
