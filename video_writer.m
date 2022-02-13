filepath = "Earth_nadir_600us/";
fps = 1;
refine_factor = 3;
load(filepath + "Ampraw.mat")

%Ampraw = Ampraw_avg;

[pxX,pxY,nbrOfFrames] = size(Ampraw);
video=zeros(120,160,nbrOfFrames);
Vq = zeros((pxX-1)*2^refine_factor+1,(pxY-1)*2^refine_factor+1,nbrOfFrames);

%%% Videos
v = VideoWriter(char(filepath + 'Raw.avi'),'Grayscale AVI');
v.FrameRate = fps;
v2 = VideoWriter(char(filepath + 'Interpolated.avi'),'Grayscale AVI');
v2.FrameRate = fps;

for k=1:nbrOfFrames
    disp(k/nbrOfFrames*100)
    video(:,:,k) = Ampraw(:,:,k);
    Vq(:,:,k) = interp2(video(:,:,k),refine_factor,'makima');
end

video = video - min(min(min(video)));
video = video/(max(max(max(video))));
Vq = Vq - min(min(min(Vq)));
Vq = Vq/(max(max(max(Vq))));
    
open(v)
writeVideo(v,min(video,1))
close(v)
open(v2)
writeVideo(v2,min(Vq,1))
close(v2)


figure
imagesc(video(:,:,1,1));
colormap gray
axis image
axis off
title('Original Image');

figure
imagesc(Vq(:,:,1));
colormap gray
axis image
axis off
title('Interpolation');

% %%% Video 2
% v2 = VideoWriter(char(filepath + 'video_10fps.avi'),'Grayscale AVI');
% v2.FrameRate = 10;
% 
% open(v2)
% for k=1:3:nbrOfFrames
%    writeVideo(v2,min(video(:,:,1,k),1))
% end
% close(v2)
% 
% 
% %%% Video 3
% v3 = VideoWriter(char(filepath + 'video_20fps.avi'),'Grayscale AVI');
% v3.FrameRate = 20;
% 
% open(v3)
% for k=1:3:nbrOfFrames
%    writeVideo(v3,min(video(:,:,1,k),1))
%    writeVideo(v3,min(video(:,:,1,k+1),1))
% end
% close(v3)