function FramesFromVideo(InputPathAndName, OutputFolder, OutputName, OutputFormat, ImageStep)
% AUTHOR: Filippo Piccinini (E-mail: f.piccinini@unibo.it)
% DATE: 03 June 2013
% NAME: FramesFromVideo
% 
% To extract and to save the single frams of a image video. The input video
% must be recorded with no compression and no interlacing. 
%
% PARAMETERS:
%  InputPathAndName Absolute or relative path where the video is stored in
%                   the computer. The Path is follow by the name of the 
%                   video. accepted formats for the input video, are 
%                   '.avi', '.flv' and '.wmv'. E.g.: 
%                   'C:\AbsolutePath\INPUTVIDEO\VideoMesenchymal.avi'.
%  OutputFolder     Absolute or relative path where the output images must 
%                   be saved in the computer. E.g.: 'OUTPUTIMAGES\' 
%  OutputName       Base name (without the format) selected for the 
%                   output images. E.g.: 'ImageMesenchymal'.
%  OutputFormat     Format selected for the output images with the dot 
%                   before the format name. E.g.: '.tif'.
%  ImageStep        1 to read all the image of the video. 2 to read the
%                   image 1, 3, 5... . Etc.
%
% OUTPUT:
%  The frames extracted from the video are saved in the OutputFolder 
%  according to the OutputName and OutputFormat.
%
% Example of usage:
% FramesFromVideo('C:\AbsolutePath\INPUTVIDEO\VideoMesenchymal.avi', 'OUTPUTIMAGES\', 'ImageMesenchymal', 'tif', 1)

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


% check on parameter "OutputFolder":
if isunix()
    Slash = '/';
else
    Slash = '\';
end;
if ~isempty(OutputFolder)
    if ~strcmp(OutputFolder(end), Slash)
        OutputFolder = strcat(OutputFolder,Slash);
    end
end

% check on parameter "OutputFormat":
if strcmp(OutputFormat(1), '.')
    OutputFormat = OutputFormat(2:end);
end

% Every time the folder OutputFolder is deleted with its content.
[stat, mess, id] = rmdir(OutputFolder,'s');
mkdir(OutputFolder)

% Video reading.
obj = mmreader(InputPathAndName);
numFrames = get(obj, 'NumberOfFrames');
if ~isempty(obj) && isempty(numFrames)
    numFrames = 10000; % We set to "numFrames" the maximum number of frames extracted from the video. This parameter can be changed.
    disp(['It is not possible to determine the number of frames in the video. We set to 10000 the maximum number of extracted frames.'])
else
    disp(['The number of frames in the video is: ' num2str(numFrames)])
end
CharactersNumber = length(num2str(numFrames));
NumberCharactersNumber = ['%.' num2str(CharactersNumber) 'd'];

% Image writing.
j = 1;
for i = 1:ImageStep:numFrames
    disp(['Frames from video. Process frame: ' num2str(i) '/' num2str(numFrames)])
    try        
        Immagine = read(obj, i); %For reading a single frame.
    catch ME1
        break
    end
    strnum = sprintf(NumberCharactersNumber,j);
    imwrite(Immagine,[OutputFolder OutputName '_' strnum '.' OutputFormat], OutputFormat)
    j = j+1;
end