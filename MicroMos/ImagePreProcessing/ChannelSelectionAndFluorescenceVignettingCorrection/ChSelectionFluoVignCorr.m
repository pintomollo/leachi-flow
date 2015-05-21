%% Channel Selection And Fluorescence Vignetting Correction

% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 22 August 2013
% NAME: ChSelectionFluoVignCorr
% 
% To correct the vignetting effect with the possibility of selecting also 
% a single channel in case of input RGB images.

% CVG (Computer Vision Group) Toolbox
% Copyright © 2012 Filippo Piccinini, Alessandro Bevilacqua, 
% Advanced Research Center on Electronic Systems (ARCES), 
% University of Bologna, Italy. All rights reserved.
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License version 2 (or higher) 
% as published by the Free Software Foundation. This program is 
% distributed WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
% General Public License for more details.

%% INPUT PARAMETERS

ImageFolder = 'C:\...\INPUTIMAGES\'; % absolute or relative path where the images are stored inside the computer. E.g.: 'INPUTIMAGES\' or 'C:\AbsolutePath\INPUTIMAGES\'.  
ImageBaseName = 'NameOfTheImages_'; % name of the images to be processed, without the final cardinal number. The last character must be the underscore. E.g.: 'ImageMesenchymal_' for the image: "ImageMesenchymal_001.tif". The images must be saved (or copied) in the following format: "Name_#...##.format".
ChannelToBeAnalyzed = 1; % accepted values: 1 = RED; 2 = GREEN; 3 = BLU; 0 = GRAY. If the input images are RGB, only the channel selected will be processed and saved in output.
BGImatrix_folderAndName = 'C:\...\INPUTIMAGES\BGImatrix.mat'; % absolute or relative path where the Matlab matrix of the BackGround referred to the input image (e.g., "BGImatrix.mat") is stored inside the computer. The field can be also empty (in this case write: BGImatrix_folderAndName = []). 
VFmatrix_folderAndName = 'C:\...\INPUTIMAGES\VFmatrix.mat'; % absolute or relative path where the Matlab matrix of the Vignetting Function (e.g., "VFmatrix.mat") is stored inside the computer. The field can be also empty (in this case write: VFmatrix_folderAndName = []).
BGVmatrix_folderAndName = []; % absolute or relative path where the Matlab matrix of the BackGround referred to the Vignetting function is stored inside the computer. The field can be also empty (in this case write: BGVmatrix_folderAndName = []). 


%% INTERNAL PARAMETERS

NumberOfImages = []; % [] means all the image in the input folder. Othervise the number of images to be processed must mbe specified.
OutputFolder = 'ImagesProcessed\';


%% CHECKS ON THE PARAMETERS: ABSOLUTELY DO NOT CHANGE FROM HERE.

% check on "ChannelToBeAnalyzed":
if ChannelToBeAnalyzed~=0 && ChannelToBeAnalyzed~=1 && ChannelToBeAnalyzed~=2 && ChannelToBeAnalyzed~=3
    error('Invalid input value for the parameter: ChannelToBeAnalyzed.');
end

% check on "ImageBaseName":
if ~strcmp(ImageBaseName(end), '_')
    ImageBaseName = strcat(ImageBaseName,'_');
end
if ~strcmp(ImageBaseName(end), '*')
    ImageBaseName = strcat(ImageBaseName,'*');
end

% check on "ImageFolder":
if isunix()
    Slash = '/';
else
    Slash = '\';
end;
if ~strcmp(ImageFolder(end), Slash)
    ImageFolder = strcat(ImageFolder,Slash);
end
if ~strcmp(OutputFolder(end), Slash)
    OutputFolder = strcat(OutputFolder,Slash);
end


%% IMAGE PROCESSING

% Every time the folder OutputFolder is deleted with its content.
[stat, mess, id] = rmdir(OutputFolder,'s');
mkdir(OutputFolder)

% Load BackGround matrix referred to the input image: the background matrix must be of the same size of the input images.
BGIfield = [];
if ~isempty(BGImatrix_folderAndName)
    BGIfieldStructure = load(BGImatrix_folderAndName);
    BGIfieldCell = struct2cell(BGIfieldStructure);
    for i = 1:length(BGIfieldCell)
        if isnumeric(BGIfieldCell{i})
            BGIfield = double(BGIfieldCell{i});
            break
        end
    end
    clear BGIfieldStructure BGIfieldCell
end

% Load BackGround matrix referred to the vignetting function: the background matrix must be of the same size of the input images.
BGVfield = [];
if ~isempty(BGVmatrix_folderAndName)
    BGVfieldStructure = load(BGVmatrix_folderAndName);
    BGVfieldCell = struct2cell(BGVfieldStructure);
    for i = 1:length(BGVfieldCell)
        if isnumeric(BGVfieldCell{i})
            BGVfield = double(BGVfieldCell{i});
            break
        end
    end
    clear BGVfieldStructure BGVfieldCell
end

% Load Vignetting Function matrix: the background matrix must be of the same size of the input images.
VGfield = [];
if ~isempty(VFmatrix_folderAndName)
    VFfieldStructure = load(VFmatrix_folderAndName);
    VFfieldCell = struct2cell(VFfieldStructure);
    for i = 1:length(VFfieldCell)
        if isnumeric(VFfieldCell{i})
            VFfield = double(VFfieldCell{i});
            break
        end
    end
    clear VFfieldStructure VFfieldCell
end

% Setting of the internal parameter based on the input images.
dirList = dir([ImageFolder ImageBaseName]);
dirList_length = length(dirList);
if isempty(NumberOfImages)
    NumberOfImages = dirList_length;
end
if dirList_length < NumberOfImages
    NumberOfImages = dirList_length;
end
disp(['Number of processed images: ' num2str(NumberOfImages)])

StringNameImage1 = dirList(1,1).name;
PositionsUnderScores = strfind(StringNameImage1, '_');
PositionLastUnderScores = PositionsUnderScores(end);
PositionsPoints = strfind(StringNameImage1, '.');
PositionLastPoint = PositionsPoints(end);
ImageFormat = StringNameImage1(PositionLastPoint:end);
ImageNameWoFormat = StringNameImage1(1:PositionLastPoint-1);
CharactersNumber = StringNameImage1(PositionLastUnderScores+1:PositionLastPoint-1);
NumberCharactersNumber = ['%.' num2str(length(CharactersNumber)) 'd'];

% Loading of the first image.
First = imread([ImageFolder dirList(1).name]);
[row, col, ch] = size(First);
if isa(First, 'uint16'); First  = 255.*double(First)./(2^16-1); end
if isa(First, 'uint32'); First  = 255.*double(First)./(2^32-1); end
if isa(First, 'uint64'); First  = 255.*double(First)./(2^64-1); end
if ch == 1
    if ChannelToBeAnalyzed == 1 || ChannelToBeAnalyzed == 2 || ChannelToBeAnalyzed == 3
        disp('The input: ChannelToBeAnalyzed has been not considered because the input images are monochannel.')
    end
    ChannelToBeAnalyzed = 1;
    ChannelAnalysed = 1;
end
if ChannelToBeAnalyzed == 1
    ChannelAnalysed = 1;
elseif ChannelToBeAnalyzed == 2
    ChannelAnalysed = 2;
elseif ChannelToBeAnalyzed == 3
    ChannelAnalysed = 3;
elseif ChannelToBeAnalyzed == 0
    ChannelAnalysed = 1;
end

% Error check.
if isempty(BGIfield)
    BGIfield = double(zeros(row, col));
end
if isempty(BGVfield)
    BGVfield = double(zeros(row, col));
end
if isempty(VFfield)
    VFfield = double(ones(row, col));
end
if size(BGIfield)~=size(BGVfield)
    error('The BackGround matrices must be of the same size.')
end
if size(BGIfield)~=size(VFfield)
    error('The BackGround matrix and the Vignetting Function matrix must be of the same size.')
end
if (size(BGIfield,1)~=size(First,1)) || (size(BGIfield,2)~=size(First,2))
    error('The BackGround matrix and the Vignetting Function matrix must be of the same size of the input images.')
end

% Channel selection and vignetting correction
for k = 1:NumberOfImages
    % Loading of the image k.
    Input = imread([ImageFolder dirList(k).name]);
    
    if isa(Input, 'uint16'); Input  = 255.*double(Input)./(2^16-1); end
    if isa(Input, 'uint32'); Input  = 255.*double(Input)./(2^32-1); end
    if isa(Input, 'uint64'); Input  = 255.*double(Input)./(2^64-1); end
    
    % Channel selection
    IndNaN = find(isnan(Input(:,:,1)));
    if ChannelToBeAnalyzed == 0
        Input = rgb2gray(Input);
    end
    Original = Input(:,:,ChannelAnalysed);
    Original = double(Original);
    NormFact = VFfield - BGVfield;
    
    % Vignetting correction
    ImgPreProcessed = ((Original - BGIfield)./NormFact).*mean(NormFact(:)) + mean(BGIfield(:));
    Ind255 = find(ImgPreProcessed>255);
    Ind0 = find(ImgPreProcessed<0);
    ImgPreProcessed(Ind255) = 255;
    ImgPreProcessed(Ind0) = 0;
    
    % Output saving
    disp(['Image number: ' num2str(k) '/' num2str(NumberOfImages) '. Channel analised: ' num2str(ChannelToBeAnalyzed)])
    StringNameImagek = dirList(k).name;
    ImageNameWoFormat = StringNameImagek(1:PositionLastPoint-1);
    Name = [OutputFolder ImageNameWoFormat '.mat'];
    save(Name, 'ImgPreProcessed')
    
    clear Input Corrected Original;
end

%% OUTPUT SAVING
copyfile('ChSelectionFluoVignCorr.m', [OutputFolder 'ChSelectionFluoVignCorr.m']);