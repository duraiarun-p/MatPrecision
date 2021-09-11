function [geomtform,Rfixed,Rmoving,optimizer,metric]=dblmrigidim(fixedVolume,movingVolume,fixedHeader,movingHeader)
[optimizer,metric] = imregconfig('multimodal');

Rfixed  = imref3d(size(fixedVolume),fixedHeader.PixelSpacing(2),fixedHeader.PixelSpacing(1),fixedHeader.SliceThickness);
Rmoving = imref3d(size(movingVolume),movingHeader.PixelSpacing(2),movingHeader.PixelSpacing(1),movingHeader.SliceThickness);

optimizer.InitialRadius = 0.004;
optimizer.MaximumIterations = 200;
geomtform = imregtform(movingVolume,Rmoving, fixedVolume,Rfixed, 'rigid', optimizer, metric);
end