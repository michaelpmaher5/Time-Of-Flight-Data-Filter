% Temporary script

Norm = zeros(120,160,nbrOfFrames);
for k=1:nbrOfFrames
    for j=1:160
        for i=1:120
            Norm(i,j,k) = norm([X(i,j,k) Y(i,j,k) Z(i,j,k)]);
        end
    end
end

Mean_Norm = zeros(120,160);
for j=1:160
    for i=1:120
        Mean_Norm(i,j) = mean(Norm(i,j,find(Norm(i,j,:))));
    end
end

Std_Norm = zeros(120,160);
for j=1:160
    for i=1:120
        Std_Norm(i,j) = std(Norm(i,j,find(Norm(i,j,:))),1);
    end
end

Mean_big = zeros(120,160,3);
for j=1:160
    for i=1:120
        Mean_big(i,j,1) = mean(X(i,j,find(Z(i,j,:))));
        Mean_big(i,j,2) = mean(Y(i,j,find(Z(i,j,:))));
        Mean_big(i,j,3) = mean(Z(i,j,find(Z(i,j,:))));
    end
end

Norm_Mean = zeros(120,160);
for j=1:160
    for i=1:120
        Norm_Mean(i,j) = norm([Mean_big(i,j,1) Mean_big(i,j,2) Mean_big(i,j,3)]);
    end
end

Std_big = zeros(120,160,3);
for j=1:160
    for i=1:120
        Std_big(i,j,1) = std(X(i,j,find(Z(i,j,:))),1);
        Std_big(i,j,2) = std(Y(i,j,find(Z(i,j,:))),1);
        Std_big(i,j,3) = std(Z(i,j,find(Z(i,j,:))),1);
    end
end

Norm_Std = zeros(120,160);
for j=1:160
    for i=1:120
        Norm_Std(i,j) = norm([Std_big(i,j,1) Std_big(i,j,2) Std_big(i,j,3)]);
    end
end

ErrMean = abs(Mean_Norm - Norm_Mean);
ErrStd = abs(Std_Norm - Norm_Std);

%%

Norm2 = zeros(13,13,nbrOfFrames);
for k=1:nbrOfFrames
    for j=1:13
        for i=1:13
            Norm2(i,j,k) = norm([M_mean_pos(i,j,1,k) M_mean_pos(i,j,2,k) M_mean_pos(i,j,3,k)]);
        end
    end
end

Mean_Norm2 = zeros(13,13);
for j=1:13
    for i=1:13
        Mean_Norm2(i,j) = mean(Norm2(i,j,:));
    end
end

Std_Norm2 = zeros(13,13);
for j=1:13
    for i=1:13
        Std_Norm2(i,j) = std(Norm2(i,j,:),1);
    end
end

Norm_Mean2 = zeros(13,13);
for j=1:13
    for i=1:13
        Norm_Mean2(i,j) = norm([Mean(i,j,1) Mean(i,j,2) Mean(i,j,3)]);
    end
end

Norm_Std2 = zeros(13,13);
for j=1:13
    for i=1:13
        Norm_Std2(i,j) = norm([Std(i,j,1) Std(i,j,2) Std(i,j,3)]);
    end
end

ErrMean2 = abs(Mean_Norm2 - Norm_Mean2);
ErrStd2 = abs(Std_Norm2 - Norm_Std2);

%%

Delta = zeros(120,160);
for j=1:160
    for i=1:120
        if max(Z(i,j,:))~=0 && min(Z(i,j,:))~=0 
            Delta(i,j) = max(Norm(i,j,:))-min(Norm(i,j,:));
        end
    end
end
 
Delta2 = zeros(13,13);
for j=1:13
    for i=1:13
        Delta2(i,j) = max(Norm2(i,j,:))-min(Norm2(i,j,:));
    end
end

%%

figure('Units','normalized','Position',[0 0 1 1]);
hold on
grid on
daspect([1 1 1])
for k=1:nbrOfFrames
    plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos(:,:,k),'.');
end

figure('Units','normalized','Position',[0 0 1 1]);
hold on
grid on
daspect([1 1 1])
for k=1:nbrOfFrames
    plot3(M_mean_pos(:,:,1,k),-M_mean_pos(:,:,2,k),M_mean_pos(:,:,3,k),'.');
end

figure('Units','normalized','Position',[0 0 1 1]);
hold on
grid on
daspect([1 1 1])
for k=1:nbrOfFrames
    plot3(Xpos(:,:,k),-Ypos(:,:,k),Norm(:,:,k),'.');
end

figure('Units','normalized','Position',[0 0 1 1]);
hold on
grid on
daspect([1 1 1])
for k=1:nbrOfFrames
    plot3(M_mean_pos(:,:,1,k),-M_mean_pos(:,:,2,k),Norm2(:,:,k),'.');
end




%%
        
figure('Units','normalized','Position',[0 0 1 1]);
plot3(Xpos(:,:,1),-Ypos(:,:,1),Delta,'.');
hold on
surf(Xpos(:,:,1),-Ypos(:,:,1),Delta,'FaceColor','interp')
daspect([1 1 1])
grid on

figure('Units','normalized','Position',[0 0 1 1]);
plot3(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Delta2,'.');
hold on
surf(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Delta2,'FaceColor','interp')
daspect([1 1 1])
grid on

%%

figure('Units','normalized','Position',[0 0 1 1]);
plot3(Xpos(:,:,1),-Ypos(:,:,1),Mean_Norm,'.');
hold on
surf(Xpos(:,:,1),-Ypos(:,:,1),Mean_Norm,'FaceColor','interp')
%daspect([1 1 1])
grid on
% 
% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(Xpos(:,:,1),-Ypos(:,:,1),Norm_Mean,'.');
% hold on
% surf(Xpos(:,:,1),-Ypos(:,:,1),Norm_Mean,'FaceColor','interp')
% %daspect([1 1 1])
% grid on

% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(Xpos(:,:,1),-Ypos(:,:,1),ErrMean,'.');
% hold on
% surf(Xpos(:,:,1),-Ypos(:,:,1),ErrMean,'FaceColor','interp')
% %daspect([1 1 0.000001])
% grid on


figure('Units','normalized','Position',[0 0 1 1]);
plot3(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Mean_Norm2,'.');
hold on
surf(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Mean_Norm2,'FaceColor','interp')
%daspect([1 1 0.5])
grid on
% 
% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Norm_Mean2,'.');
% hold on
% surf(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Norm_Mean2,'FaceColor','interp')
% %daspect([1 1 0.5])
% grid on

% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),ErrMean2,'.');
% hold on
% surf(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),ErrMean2,'FaceColor','interp')
% %daspect([1 1 0.001])
% grid on

%%

figure('Units','normalized','Position',[0 0 1 1]);
plot3(Xpos(:,:,1),-Ypos(:,:,1),Std_Norm,'.');
hold on
surf(Xpos(:,:,1),-Ypos(:,:,1),Std_Norm,'FaceColor','interp')
%daspect([1 1 1])
grid on
% 
% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(Xpos(:,:,1),-Ypos(:,:,1),Norm_Std,'.');
% hold on
% surf(Xpos(:,:,1),-Ypos(:,:,1),Norm_Std,'FaceColor','interp')
% %daspect([1 1 1])
% grid on

% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(Xpos(:,:,1),-Ypos(:,:,1),ErrStd,'.');
% hold on
% surf(Xpos(:,:,1),-Ypos(:,:,1),ErrStd,'FaceColor','interp')
% %daspect([1 1 0.000001])
% grid on


figure('Units','normalized','Position',[0 0 1 1]);
plot3(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Std_Norm2,'.');
hold on
surf(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Std_Norm2,'FaceColor','interp')
%daspect([1 1 0.5])
grid on
% 
% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Norm_Std2,'.');
% hold on
% surf(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),Norm_Std2,'FaceColor','interp')
% %daspect([1 1 0.5])
% grid on

% figure('Units','normalized','Position',[0 0 1 1]);
% plot3(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),ErrStd2,'.');
% hold on
% surf(M_mean_pos(:,:,1,1),-M_mean_pos(:,:,2,1),ErrStd2,'FaceColor','interp')
% %daspect([1 1 0.001])
% grid on

