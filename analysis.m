clear, clc, close all
%% Part 1
dicomImage = dicomread("image10.dcm");
figure
imshow(dicomImage, [])
title('Original DICOM Image')

croppedImage = imcrop(dicomImage,[175 125 127 127]); % Crop 128x128 section
figure
imshow(croppedImage, [])
title('Cropped 128x128 Section')
%% Part 2
numAngles = 180; 
projections = zeros(numAngles, size(croppedImage, 1));

for angle = 1:180
    rotatedImage = imrotate(croppedImage, angle, 'bilinear', 'crop'); 
    projections(angle, :) = sum(rotatedImage, 1); 
end

figure
imagesc(projections)
colormap('gray')
xlabel('Detector Pixels')
ylabel('Projection Angle (degrees)')
title('180 1D Projections (Sinogram)')
%% Part 3
pad=padarray(croppedImage,[50 50],'both');
for i = 1:180
    rotate=imrotate(pad,(i-1),'bilinear','crop');
    sinogram(:,i) = sum(rotate);
end
dfft = fftshift(fft(fftshift(sinogram)));
figure
imshow(log(abs(dfft)),[])
%% Part 4 and 5
theta = repmat((0:1:(180-1))*pi/180,length(dfft),1);
bound = (length(dfft)-1)/2;
rho = repmat((-ceil(bound):floor(bound))',1,180);
[x,y] = pol2cart(theta,rho);
[XI,YI] = meshgrid((-ceil(bound):floor(bound)));
assembled_dfft = griddata(x,y,dfft,XI,YI,'linear');
assembled_dfft(isnan(assembled_dfft))=0;

figure
imshow(log(abs(assembled_dfft)),[])
colormap('gray')
title('Interpolated 2D Fourier Transform')
%% Part 6
direct = fft2(croppedImage);
figure
imshow(log(abs(fftshift(direct))),[])
colormap('gray')
title('2D Fourier Transform (Magnitude)')
%% Part 8
inverseImage = abs(fftshift(ifft2(fftshift(assembled_dfft)))');
figure
imagesc(inverseImage)
colormap('gray')
title('Reconstructed Image from Inverse 2D Fourier Transform')
