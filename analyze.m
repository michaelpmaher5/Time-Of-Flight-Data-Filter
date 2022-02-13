close all

filepath = "200_5/";                                %% File data is being pulled from
fps = 5;                                            %% FPS

X_refine_index = 1+12*4;                            %% X length of each square on sim
Y_refine_index = 1+12*4;

load(filepath + "Xraw.mat");                         %% XPos of sail at each point
load(filepath + "Yraw.mat");                         %% YPos of sail at each point
load(filepath + "Zraw.mat");                         %% ZPos of sail at each point
load(filepath + "Ampraw.mat");                       %% Displacement of sail at each point wrt reference

[pxX,pxY,nbrOfFrames] = size(Ampraw);
%nbrOfFrames=2;

%%% We get rid of the values where Z=-1 %%%
Xpos = Xraw;                X = Xraw;
Ypos = Yraw;                Y = Yraw;
Zpos = Zraw;                Z = Zraw;
Amppos = Ampraw;            Amp = Ampraw;

index = (Zraw==-1);
Xpos(index) = nan;          X(index) = 0;    %%pos is used for plotting, nan to make it clean
Ypos(index) = nan;          Y(index) = 0;    %%X,Y,Z is used for calculation
Zpos(index) = nan;          Z(index) = 0;
%Amppos(index) = nan;        Amp(index) = 0;
index = (Zraw==0);
Xpos(index) = nan;          X(index) = 0;
Ypos(index) = nan;          Y(index) = 0;
Zpos(index) = nan;          Z(index) = 0;
%Amppos(index) = nan;        Amp(index) = 0;


%------------------------------------------------------------------------------------------------------
%% Analysis %%%

% Meshgrid (Plotting only)
ix = 1:160; %% Length of X side of mesh
iy = 1:120; %% Length of Y side of mesh
[iX,iY] = meshgrid(ix,iy);
[Xp,Yp] = meshgrid(30:10:130,5:10:105); 


% We get the outer corners and the sides' lengths
[S1,S2,S3,S4] = summit(Z,nbrOfFrames);  %%detect the location of each corner in pixels (3,1,62) (xyz at each frame)

[ a,b,c,d ] = plane( S1,S2,S3,S4,nbrOfFrames );  %% a b c d is the coefficients for average plane equation

[ index,limit ] = value_limit( Z,a,b,c,d,nbrOfFrames );  %% removing points outside of plane boundary
Xpos(index) = nan;          X(index) = 0;
Ypos(index) = nan;          Y(index) = 0;
Zpos(index) = nan;          Z(index) = 0;
Amppos(index) = nan;        Amp(index) = 0;

% [Z_flat] = flattener(Z,a,b,c,d,nbrOfFrames);
% [Zpos_flat] = flattener(Zpos,a,b,c,d,nbrOfFrames);
[Z_flat] = flattener_v2(Z,a,b,c,d,nbrOfFrames);             %% since the camera was tilted, this is to rotate the model to look level
[Zpos_flat] = flattener_v2(Zpos,a,b,c,d,nbrOfFrames);       %% returns plane completely rotated

[ M_edges,E1,E2,E3,E4 ] = edges( S1,S2,S3,S4,Z,nbrOfFrames );   %% The red dots on graph, stored as x y z pos for each frame along respective axis
M_edges(M_edges==0)=NaN;                %% removed 0 to make plot look nice

[ L1,L2,L3,L4 ] = lengths( E1,E2,E3,E4,nbrOfFrames );    %% length of each side of sail at each frame
 
[ M_bound, A1, A2, A3, A4 ] = boundery_matrix( S1, S2, S3, S4, nbrOfFrames );  %% not used 
[ M_bound_v2, E1_bound, E2_bound, E3_bound, E4_bound ] = boundery_matrix_v2( L1,L2,L3,L4,E1,E2,E3,E4,nbrOfFrames ); %% Mboundv2 is only used, returns coordinates for the red crosses for the grid

[ M_mean_pos ] = mean_matrix_v2( M_bound_v2, X, Y, Z, nbrOfFrames );    %% 4-D: 1st and 2nd are indexes (x y), 3rd is 3 values ( pos x y and z), 4th is frames
[ M_mean_pos_flat ] = flattener_v2(M_mean_pos,a,b,c,d,nbrOfFrames);     %% same flattener thing to make it level

Mean = mean_global_matrix(M_mean_pos,nbrOfFrames);     %% Average of the drag sail amplitudes at each pos over each frame
Std = std_matrix(M_mean_pos,nbrOfFrames);              %% std dev over each frame
Mean_flat = mean_global_matrix(M_mean_pos_flat,nbrOfFrames); %% same as first two but using level plane
Std_flat = std_matrix(M_mean_pos_flat,nbrOfFrames);

refV = getVelAcc(M_mean_pos_flat,fps,nbrOfFrames);
refA = getVelAcc(refV,fps,nbrOfFrames);

refVstddev = standarddeviation(refV,nbrOfFrames);
refAstddev = standarddeviation(refA,nbrOfFrames);

%% Kalman Filter
stddev = standarddeviation(M_mean_pos_flat,nbrOfFrames);
filteredpos = kalmanfilter(M_mean_pos_flat,stddev,nbrOfFrames);
filteredposstddev = standarddeviation(filteredpos,nbrOfFrames);
filteredV = getVelAcc(filteredpos,fps,nbrOfFrames);
estimatedA = getVelAcc(filteredV,fps,nbrOfFrames);

Vstddev = standarddeviation(filteredV,nbrOfFrames);
kalmanV = kalmanfilter(filteredV,Vstddev,nbrOfFrames);

filteredA = getVelAcc(kalmanV,fps,nbrOfFrames);
Astddev = standarddeviation(filteredA,nbrOfFrames);
kalmanA = kalmanfilter(filteredA,Astddev,nbrOfFrames);



%% Interpolation
Mq = interpolation(M_mean_pos,X_refine_index,Y_refine_index,nbrOfFrames);   
%Meanq = interpolation(Mean,X_refine_index,Y_refine_index,1);
Meanq = interpolation(Mean_flat,X_refine_index,Y_refine_index,1);

% Raw Interpolation
[Xrawposmean,Yrawposmean,Zrawposmean] = rawinterpolate(Xraw,Yraw,Zraw,nbrOfFrames);

% Time derivatives
[ V,A ] = time_derivatives( Mq,fps,nbrOfFrames );
V_mean = mean_global_matrix(V,nbrOfFrames-1);
A_mean = mean_global_matrix(A,nbrOfFrames-2);



M_mean_pos(7,7,:,:) = nan;   %% makes the white square in the middle because we dont care about it
Mean(7,7,:) = nan;
M_mean_pos_flat(7,7,:,:) = nan;
Mean_flat(7,7,:) = nan;
%% Surface Norm
normalsfiltered = getNormalCoordinates(filteredpos,nbrOfFrames);
%% Angular Values


polar = converttopolar(normalsfiltered,nbrOfFrames);
polarV = getVelAcc(polar,fps,nbrOfFrames);
polarA = getVelAcc(polarV,fps,nbrOfFrames);

originalpolar = converttopolar(M_mean_pos_flat,nbrOfFrames);
originalpolarV = getVelAcc(originalpolar,fps,nbrOfFrames);
originalpolarA = getVelAcc(originalpolarV,fps,nbrOfFrames);

originalpolarstddev = standarddeviation(originalpolar,nbrOfFrames);
originalpolarVstddev = standarddeviation(originalpolarV,nbrOfFrames);
originalpolarAstddev = standarddeviation(originalpolarA,nbrOfFrames);

polarVstddev = standarddeviation(polarV,nbrOfFrames);
kalmanpolarV = kalmanfilter(polarV,polarVstddev,nbrOfFrames);
kalmanpolarA = getVelAcc(kalmanpolarV,fps,nbrOfFrames);

stddevPolarPFiltered = standarddeviation(polar,nbrOfFrames);
stddevPolarAFiltered = standarddeviation(kalmanpolarA,nbrOfFrames);

%% Average Standard Deviations
aveStDevPosUnfiltered = getAverageStandardDeviation(stddev); %unfiltered values
aveStDevVelUnfiltered = getAverageStandardDeviation(refVstddev);
aveStDevAccUnfiltered = getAverageStandardDeviation(refAstddev);

aveStDevP = getAverageStandardDeviation(filteredposstddev); %filtered values
aveStDevV = getAverageStandardDeviation(Vstddev);
aveStDevA = getAverageStandardDeviation(Astddev);

aveStDevPolarPosUnfiltered = getAverageStandardDeviation(originalpolarstddev);
aveStDevPolarVelUnfiltered = getAverageStandardDeviation(originalpolarVstddev);
aveStDevPolarAccUnfiltered = getAverageStandardDeviation(originalpolarAstddev);

aveStDevPolarPos = getAverageStandardDeviation(stddevPolarPFiltered);
aveStDevPolarVel = getAverageStandardDeviation(polarVstddev);
aveStDevPolarAcc = getAverageStandardDeviation(stddevPolarAFiltered);

values = [aveStDevPosUnfiltered aveStDevP; aveStDevVelUnfiltered aveStDevV; aveStDevAccUnfiltered aveStDevA; aveStDevPolarPosUnfiltered aveStDevPolarPos; aveStDevPolarVelUnfiltered aveStDevPolarVel; aveStDevPolarAccUnfiltered aveStDevPolarAcc];
%------------------------------------------------------------------------------------------------------
%% Plot %%%


f0=figure('Units','normalized','Position',[0 0 1 1],'Name','Plot point');  

hold on;
for k=1:1
    clf
    
    subplot(2,2,1)
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
     hold on;
    surf(M_mean_pos(:,:,1,k),-M_mean_pos(:,:,2,k),M_mean_pos(:,:,3,k),'FaceColor','interp');
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    colorbar()
    daspect([1 1 1])
    xlabel('X Position (mm)')
ylabel('Y Position (mm)')
zlabel('Z Position (mm)')
    grid on
    
    subplot(2,2,3)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    surf(M_mean_pos_flat(:,:,1,k),-M_mean_pos_flat(:,:,2,k),M_mean_pos_flat(:,:,3,k),'FaceColor','interp');
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
        xlabel('X Position (mm)')
ylabel('Y Position (mm)')
zlabel('Z Position (mm)')
    colorbar()
    daspect([1 1 1])
    grid on
    
    
    
    subplot(1,2,2)
    hold on;
    plot3(iX,-iY,Zpos(:,:,k),'.b');
    daspect([1 1 50])
    xlim([30 135]);
    ylim([-110 -5]);
	zlim([-10 4000]);
        xlabel('X Position (mm)')
ylabel('Y Position (mm)')
zlabel('Z Position (mm)')
    grid on
    
    %Offset planes
%     line(Xp,-Yp,-1/c(k)*(a(k)*Xp+b(k)*Yp+d(k))+limit(k),'Color','y');
%     line(Xp,-Yp,-1/c(k)*(a(k)*Xp+b(k)*Yp+d(k))-limit(k),'Color','y');
%     line(Xp.',-Yp.',-1/c(k)*(a(k)*Xp+b(k)*Yp+d(k)).'+limit(k),'Color','y');
%     line(Xp.',-Yp.',-1/c(k)*(a(k)*Xp+b(k)*Yp+d(k)).'-limit(k),'Color','y');
    
    %Summits
    plot3(S1(1,1,k),-S1(2,1,k),S1(3,1,k),'or');
    plot3(S2(1,1,k),-S2(2,1,k),S2(3,1,k),'or');
    plot3(S3(1,1,k),-S3(2,1,k),S3(3,1,k),'or');
    plot3(S4(1,1,k),-S4(2,1,k),S4(3,1,k),'or');
    
    plot3(iX,-iY,M_edges(:,:,k),'*r');
    
    line(M_bound_v2(:,:,1,k),-M_bound_v2(:,:,2,k),'Color','g'); 
    line(M_bound_v2(:,:,1,k).',-M_bound_v2(:,:,2,k).','Color','g'); 
%     line([M_bound_v2(1,:,1,k);M_bound_v2(14,:,1,k)],[-M_bound_v2(1,:,2,k);-M_bound_v2(14,:,2,k)],'Color','g'); 
% 	line([M_bound_v2(:,1,1,k).';M_bound_v2(:,14,1,k).'],[-M_bound_v2(:,1,2,k).';-M_bound_v2(:,14,2,k).'],'Color','g'); 
	
    plot(M_bound_v2(:,:,1,k),-M_bound_v2(:,:,2,k),'+r');
    
         
    pause(1/fps)
    
end
% 
% 
% 
% figure('Units','normalized','Position',[0 0 1 1],'Name','Plot point');  
% 
% hold on;
% for k=1:1
%     clf
%     
%     subplot(1,2,1)
%     hold on;
%     plot3(iX,-iY,Zpos(:,:,k),'.b');
%     daspect([1 1 100])
%     xlim([30 135]);
%     ylim([-110 -5]);
% 	zlim([-10 4000]);
%     grid on
%     
%     pS1 = plot3(S1(1,1,k),-S1(2,1,k),S1(3,1,k),'or');
%     pS2 = plot3(S2(1,1,k),-S2(2,1,k),S2(3,1,k),'or');
%     pS3 = plot3(S3(1,1,k),-S3(2,1,k),S3(3,1,k),'or');
%     pS4 = plot3(S4(1,1,k),-S4(2,1,k),S4(3,1,k),'or');
%     
%     pL1 = line([S1(1,1,k) S2(1,1,k)],[-S1(2,1,k) -S2(2,1,k)],[S1(3,1,k) S2(3,1,k)],'Color','r');
%     pL2 = line([S2(1,1,k) S3(1,1,k)],[-S2(2,1,k) -S3(2,1,k)],[S2(3,1,k) S3(3,1,k)],'Color','r');
%     pL3 = line([S3(1,1,k) S4(1,1,k)],[-S3(2,1,k) -S4(2,1,k)],[S3(3,1,k) S4(3,1,k)],'Color','r');
%     pL4 = line([S4(1,1,k) S1(1,1,k)],[-S4(2,1,k) -S1(2,1,k)],[S4(3,1,k) S1(3,1,k)],'Color','r');
%     
%     for p=1:14
%         line([A1(1,p,k) A3(1,p,k)],-[A1(2,p,k) A3(2,p,k)],[0 0],'Color','g');
%         line([A4(1,p,k) A2(1,p,k)],-[A4(2,p,k) A2(2,p,k)],[0 0],'Color','g');
%     end
%     
%     for i=1:14
%         for j=1:14
%             plot(M_bound(i,j,1,k),-M_bound(i,j,2,k),'+r');
%         end
%     end
%     
%     
%     subplot(1,2,2)
%     hold on;
%     plot3(iX,-iY,Zpos(:,:,k),'.b');
%     daspect([1 1 100])
%     xlim([30 135]);
%     ylim([-110 -5]);
% 	zlim([-10 4000]);
%     grid on
%     
%     %Summits
%     plot3(S1(1,1,k),-S1(2,1,k),S1(3,1,k),'or');
%     plot3(S2(1,1,k),-S2(2,1,k),S2(3,1,k),'or');
%     plot3(S3(1,1,k),-S3(2,1,k),S3(3,1,k),'or');
%     plot3(S4(1,1,k),-S4(2,1,k),S4(3,1,k),'or');
%     
%     plot3(iX,-iY,M_edges(:,:,k),'*r');
%     
%     line(M_bound_v2(:,:,1,k),-M_bound_v2(:,:,2,k),'Color','g'); 
%     line(M_bound_v2(:,:,1,k).',-M_bound_v2(:,:,2,k).','Color','g'); 
% %     line([M_bound_v2(1,:,1,k);M_bound_v2(14,:,1,k)],[-M_bound_v2(1,:,2,k);-M_bound_v2(14,:,2,k)],'Color','g'); 
% % 	  line([M_bound_v2(:,1,1,k).';M_bound_v2(:,14,1,k).'],[-M_bound_v2(:,1,2,k).';-M_bound_v2(:,14,2,k).'],'Color','g'); 
% 	
%     plot(M_bound_v2(:,:,1,k),-M_bound_v2(:,:,2,k),'+r');
%          
%     
%     
%     pause(1/fps)
%     
%     
% 
%  end


%{
%%% Mean
figure('Units','normalized','Position',[0 0 1 1],'Name','Average');
plot3(Mean(:,:,1),-Mean(:,:,2),Mean(:,:,3),'.b')
hold on;
surf(Mean(:,:,1),-Mean(:,:,2),Mean(:,:,3),'FaceColor','interp')
daspect([1 1 1])
grid on
colorbar()
    
% %%% STD X
% figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation X');
% plot3(Mean(:,:,1),-Mean(:,:,2),Std(:,:,1),'.b')
% hold on;
% surf(Mean(:,:,1),-Mean(:,:,2),Std(:,:,1),'FaceColor','interp')
% daspect([1 1 0.5])
% grid on
% colorbar()
% 
% %%% STD Y
% figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation Y');
% plot3(Mean(:,:,1),-Mean(:,:,2),Std(:,:,2),'.b')
% hold on;
% surf(Mean(:,:,1),-Mean(:,:,2),Std(:,:,2),'FaceColor','interp')
% daspect([1 1 0.5])
% grid on
% colorbar()
% 
% %%% STD Z
% figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation Z');
% plot3(Mean(:,:,1),-Mean(:,:,2),Std(:,:,3),'.b')
% hold on;
% surf(Mean(:,:,1),-Mean(:,:,2),Std(:,:,3),'FaceColor','interp')
% daspect([1 1 0.5])
% grid on
% colorbar()
% 
% %%% STD Total
% figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation Total');
% plot3(Mean(:,:,1),-Mean(:,:,2),Std(:,:,1)+Std(:,:,2)+Std(:,:,3),'.b')
% hold on;
% surf(Mean(:,:,1),-Mean(:,:,2),Std(:,:,1)+Std(:,:,2)+Std(:,:,3),'FaceColor','interp')
% daspect([1 1 0.5])
% grid on
% colorbar()

%%% Mean_flat
figure('Units','normalized','Position',[0 0 1 1],'Name','Average flattened');
plot3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Mean_flat(:,:,3),'.b')
hold on;
surf(Mean_flat(:,:,1),-Mean_flat(:,:,2),Mean_flat(:,:,3),'FaceColor','interp')
daspect([1 1 1])
grid on
colorbar()
    
% %%% STD_flat X
% figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation flattened X');
% plot3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,1),'.b')
% hold on;
% surf(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,1),'FaceColor','interp')
% daspect([1 1 0.5])
% grid on
% colorbar()
% 
% %%% STD_flat Y
% figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation flattened Y');
% plot3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,2),'.b')
% hold on;
% surf(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,2),'FaceColor','interp')
% daspect([1 1 0.5])
% grid on
% colorbar()
% 
%%% STD_flat Z
figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation flattened Z');
plot3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,3),'.b')
hold on;
surf(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,3),'FaceColor','interp')
daspect([1 1 0.5])
grid on
colorbar()
% 
% %%% STD_flat Total
% figure('Units','normalized','Position',[0 0 1 1],'Name','Standard Deviation flattened Total');
% plot3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,1)+Std_flat(:,:,2)+Std_flat(:,:,3),'.b')
% hold on;
% surf(Mean_flat(:,:,1),-Mean_flat(:,:,2),Std_flat(:,:,1)+Std_flat(:,:,2)+Std_flat(:,:,3),'FaceColor','interp')
% daspect([1 1 0.5])
% grid on
% colorbar()



%%% Interpolation
figure('Units','normalized','Position',[0 0 1 1],'Name','Interpolation'); 

subplot(1,2,1)
hold on;
%plot3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Mean_flat(:,:,3),'.b')
surf(Mean_flat(:,:,1),-Mean_flat(:,:,2),Mean_flat(:,:,3),'FaceColor','interp')
daspect([1 1 1])
grid on
colorbar()

subplot(1,2,2)
hold on
surf(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),'FaceColor','interp');
daspect([1 1 1])
grid on
colorbar()

    

%%% Normal vectors

% figure('Units','normalized','Position',[0 0 1 1],'Name','Normal average');
% plot3(Mean(:,:,1),-Mean(:,:,2),Mean(:,:,3),'.b')
% hold on;
% surf(Mean(:,:,1),-Mean(:,:,2),Mean(:,:,3),'FaceColor','interp');
% [U,V,W]=surfnorm(Mean(:,:,1)',-Mean(:,:,2)',Mean(:,:,3)','FaceColor','interp');
% quiver3(Mean(:,:,1),-Mean(:,:,2),Mean(:,:,3),U,V,W,1,'r')
% daspect([1 1 1])
% grid on
% colorbar()
% 
figure('Units','normalized','Position',[0 0 1 1],'Name','Normal average flattened');
plot3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Mean_flat(:,:,3),'.b')
hold on;
surf(Mean_flat(:,:,1),-Mean_flat(:,:,2),Mean_flat(:,:,3),'FaceColor','interp');
[U,V,W]=surfnorm(Mean_flat(:,:,1)',-Mean_flat(:,:,2)',Mean_flat(:,:,3)','FaceColor','interp');
quiver3(Mean_flat(:,:,1),-Mean_flat(:,:,2),Mean_flat(:,:,3),U,V,W,1,'r')
daspect([1 1 1])
grid on
colorbar()
% 
% figure('Units','normalized','Position',[0 0 1 1],'Name','Normal interpolation');
% plot3(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),'.b')
% hold on;
% surf(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),'FaceColor','interp');
% [U,V,W]=surfnorm(Meanq(:,:,1)',-Meanq(:,:,2)',Meanq(:,:,3)','FaceColor','interp');
% quiver3(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),U,V,W,1,'r')
% daspect([1 1 1])
% grid on
% colorbar()




%%% Positions

% for q=1:3:13    %column
%     
%     figure('Units','normalized','Position',[0 0 1 1],'Name',['Position, row ' int2str(q)]);  
%     
%     for p=1:3:13    %row
%         
%         subplot(5,1,(p+2)/3)
%         hold on
%         plot(squeeze(M_mean_pos(p,q,3,:)),':')
%         plot(avg_7(squeeze(M_mean_pos(p,q,3,:))),'.-')
%         %plot(sgolayfilt(squeeze(M_mean_pos(p,q,3,:)),3,7),'.-')
%         
%         grid on
%         
%     end
% end




%%% Velocities

% scale = 3;
% 
% figure('Units','normalized','Position',[0 0 1 1],'Name','Velocity');  
% 
% hold on;
% for k=1:nbrOfFrames-1
%     clf
%     
%     subplot(1,2,1)
%     hold on;
%     surf(Mq(:,:,1,k),-Mq(:,:,2,k),Mq(:,:,3,k),'FaceColor','interp');
%     daspect([1 1 1])
%     grid on
%     colorbar()
%         
%     subplot(1,2,2)
%     mesh(Mq(:,:,1,k),-Mq(:,:,2,k),Mq(:,:,3,k),'EdgeColor','interp');
%     hold on;
%     quiver3(Mq(:,:,1,k),-Mq(:,:,2,k),Mq(:,:,3,k),V(:,:,1,k),V(:,:,2,k),V(:,:,3,k),scale)
%     daspect([1 1 1])
% %     xlim([30 135]);
% %     ylim([-110 -5]);
% % 	  zlim([-10 4000]);
%     grid on
%     
%     
%     pause(1/fps)
%     
% end




% scale = 3;
% 
% figure('Units','normalized','Position',[0 0 1 1],'Name','Mean velocities and accelerations'); 
% 
% subplot(1,2,1)
% mesh(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),'EdgeColor','interp');
% hold on;
% quiver3(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),V_mean(:,:,1),V_mean(:,:,2),V_mean(:,:,3),scale)
% daspect([1 1 1])
% grid on
% 
% subplot(1,2,2)
% mesh(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),'EdgeColor','interp');
% hold on;
% quiver3(Meanq(:,:,1),-Meanq(:,:,2),Meanq(:,:,3),A_mean(:,:,1),A_mean(:,:,2),A_mean(:,:,3),scale)
% daspect([1 1 1])
% grid on









cmap=colormap(jet);
colormap(parula);
    
figure('Units','normalized','Position',[0 0 1 1],'Name','Raw values');  
hold on;
for k=1:1
    clf
    
    subplot(2,2,1)
    hold on;
    surf(iX,-iY,Xpos(:,:,k));
    daspect([1 1 100])
    colormap(cmap(end:-1:1,:))
    grid on
    
    subplot(2,2,2)
    hold on;
    surf(iX,-iY,Ypos(:,:,k));
    daspect([1 1 100])
    colormap(cmap(end:-1:1,:))
    grid on
    
    subplot(2,2,3)
    hold on;
    surf(iX,-iY,Zpos(:,:,k));
    daspect([1 1 10])
    colormap(cmap(end:-1:1,:))
    grid on
    
    subplot(2,2,4)
    hold on;
    surf(iX,-iY,Amppos(:,:,k));
    daspect([1 1 1000])
    %colormap(gray)
    grid on
    
    pause(1/fps)
end

    
f1=figure('Units','normalized','Position',[0 0 1 1],'Name','Plot point');     
plot3(Xpos(:,:,1),-Ypos(:,:,1),Zpos(:,:,1),'.b');
% xlim([-2000 2000]);
% ylim([-2000 2000]);
% zlim([2000 4000]);
daspect([1 1 1])

f2=figure('Units','normalized','Position',[0 0 1 1],'Name','Plot surface');     
surf(Xpos(:,:,1),-Ypos(:,:,1),Zpos(:,:,1));
xlim([-1500 1750]);
ylim([-1500 1750]);
% zlim([2000 4000]);
daspect([1 1 1])

f3=figure('Units','normalized','Position',[0 0 1 1],'Name','Plot point');     
plot3(Xpos(:,:,end),-Ypos(:,:,end),Zpos(:,:,end),'.b');
% xlim([-2000 2000]);
% ylim([-2000 2000]);
% zlim([2000 4000]);
daspect([1 1 1])
    
f4=figure('Units','normalized','Position',[0 0 1 1],'Name','Plot surface');     
surf(Xpos(:,:,end),-Ypos(:,:,end),Zpos(:,:,end));
%xlim([-2000 2000]);
%ylim([-1500 1500]);
%zlim([0 4000]);
daspect([1 1 1])

%}
%% Create a reference plane for graphing purposes
referenceplane = M_mean_pos_flat;

for i = 1:size(referenceplane,1)
    for j = 1:size(referenceplane,2)
        for k = 1:nbrOfFrames
            referenceplane(i,j,3,k) = 0; % Set all k values to zero so that it is just a flat plane
        end
    end
end
cmap=colormap(jet);
colormap(parula);
%% Kalman Filter Graph (Position)

figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Model Using Kalman Filter');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    surf(M_mean_pos_flat(:,:,1,k),-M_mean_pos_flat(:,:,2,k),M_mean_pos_flat(:,:,3,k),'FaceColor','interp');
    zlabel('Amplitude (mm)')
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Level Sail, Unfiltered'))
    colorbar()
    daspect([1 1 1])
    grid on
    
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    surf(filteredpos(:,:,1,k),filteredpos(:,:,2,k),filteredpos(:,:,3,k));
    zlabel('Amplitude (mm)')
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Filtered Result'))
    grid on
    colorbar()
    
    hold on;



    pause(1/fps)
end

%% Kalman Filter Graph (Velocity)

figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Velocity Model Using Kalman Filter');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
    zlabel('Amplitude (mm/s)')
    hold on;
    for x1 = 1:size(filteredpos,1)
        for y1 = 1:size(filteredpos,2)
            %put on new graph where each vector displays at relative pos
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+refV(x1,y1,1,k),filteredpos(x1,y1,2,k)+refV(x1,y1,2,k),refV(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Level Sail Velocity, Unfiltered'))
    colorbar()
    daspect([1 1 1])
    grid on
    %---------------------------------------------------------------------------------
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
    zlabel('Amplitude (mm/s)')
     hold on;


    for x1 = 1:size(filteredpos,1)
        for y1 = 1:size(filteredpos,2)
            %put on new graph where each vector displays at relative pos
            %NOTE: filteredV IS THE DERIVATIVE OF THE KALMANED POS, kalmanV
            %IS THE KALMANED filteredV!!!!!!!!!! CHANGE EACH OTHER OUT IN
            %FOLLOWING LINES TO COMPARE
            
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+filteredV(x1,y1,1,k),filteredpos(x1,y1,2,k)+filteredV(x1,y1,2,k),filteredV(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Filtered Velocity Result'))
    grid on
    colorbar()
    
    
    
    pause(1/fps)
end

%% Kalman Filter Graph (Acceleration)

figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Acceleration Model Using Kalman Filter');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
    zlabel('Amplitude (mm/s^2)')
    hold on;
    for x1 = 1:size(filteredpos,1)
        for y1 = 1:size(filteredpos,2)
            %put on new graph where each vector displays at relative pos
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+refA(x1,y1,1,k),filteredpos(x1,y1,2,k)+refA(x1,y1,2,k),refA(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
    
    title(strcat('Frame ',num2str(k),': Level Sail Acceleration, Unfiltered'))
    colorbar()
    daspect([1 1 1])
    grid on
    
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
    zlabel('Amplitude (mm/s^2)')
    hold on;
    for x1 = 1:size(filteredpos,1)
        for y1 = 1:size(filteredpos,2)
            %put on new graph where each vector displays at relative pos
            %NOTE: filteredA IS THE DERIVATIVE OF THE KALMANED VEL, kalmanA
            %IS THE KALMANED filteredA!!!!!!!!!! estimatedA ONLY HAD FILTER APPLIED AT POS MATRIX ONLY!! CHANGE EACH OTHER OUT IN
            %FOLLOWING LINES TO COMPARE
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+filteredA(x1,y1,1,k),filteredpos(x1,y1,2,k)+filteredA(x1,y1,2,k),filteredA(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Filtered Acceleration Result'))
    grid on
    colorbar()
    
    
    
    pause(1/fps)
end


%% Kalman Graphs (Angular Velocity)
figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Angular Velocity Model Using Kalman Filter');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
    zlabel('Amplitude (mm)')
    hold on;
    for x1 = 1:size(polar,1)
        for y1 = 1:size(polar,2)
            %put on new graph where each vector displays at relative pos
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+originalpolarV(x1,y1,1,k),filteredpos(x1,y1,2,k)+originalpolarV(x1,y1,2,k),originalpolarV(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Level Sail Angular Velocity, Unfiltered'))
    colorbar()
    daspect([1 1 1])
    grid on
    %---------------------------------------------------------------------------------
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
    zlabel('Amplitude (mm)')
     hold on;


    for x1 = 1:size(filteredpos,1)
        for y1 = 1:size(filteredpos,2)
            %put on new graph where each vector displays at relative pos
            %NOTE: filteredV IS THE DERIVATIVE OF THE KALMANED POS, kalmanV
            %IS THE KALMANED filteredV!!!!!!!!!! CHANGE EACH OTHER OUT IN
            %FOLLOWING LINES TO COMPARE
            
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+polarV(x1,y1,1,k),filteredpos(x1,y1,2,k)+polarV(x1,y1,2,k),polarV(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Filtered Angular Velocity Result'))
    grid on
    colorbar()
    
    
    
    pause(1/fps)
end


%% Kalman Graphs (Angular Acceleration)
figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Angular Acceleration Model Using Kalman Filter');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
    hold on;
    for x1 = 1:size(polar,1)
        for y1 = 1:size(polar,2)
            %put on new graph where each vector displays at relative pos
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+originalpolarA(x1,y1,1,k),filteredpos(x1,y1,2,k)+originalpolarA(x1,y1,2,k),originalpolarA(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Level Sail Angular Acceleration, Unfiltered'))
    colorbar()
    daspect([1 1 1])
    grid on
    %---------------------------------------------------------------------------------
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    surf(referenceplane(:,:,1,k),referenceplane(:,:,2,k),referenceplane(:,:,3,k),'FaceColor','interp');
     hold on;


    for x1 = 1:size(filteredpos,1)
        for y1 = 1:size(filteredpos,2)
            %put on new graph where each vector displays at relative pos
            %NOTE: filteredV IS THE DERIVATIVE OF THE KALMANED POS, kalmanV
            %IS THE KALMANED filteredV!!!!!!!!!! CHANGE EACH OTHER OUT IN
            %FOLLOWING LINES TO COMPARE
            
            P1 = [filteredpos(x1,y1,1,k),filteredpos(x1,y1,2,k),0]; %make z zero so that all vectors start from xy plane
            P2 = [filteredpos(x1,y1,1,k)+polarA(x1,y1,1,k),filteredpos(x1,y1,2,k)+polarA(x1,y1,2,k),polarA(x1,y1,3,k)];
            pts = [P1; P2];    
            plot3(pts(:,1), pts(:,2), pts(:,3),'-^','Color','red','MarkerSize',5)
        end
    end
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Filtered Angular Acceleration Result'))
    grid on
    colorbar()
    
    
    
    pause(1/fps)
end




%% Bar Graph Angular Velocity Original
figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Original Angular Velocity Model Bar Graph');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    b = bar3(originalpolarV(:,:,2,k));
    hold on;
    
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Angular Velocity Comp Thetadot Original'))
    colorbar()
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    daspect([1 1 1])
    grid on
    %---------------------------------------------------------------------------------
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    b = bar3(originalpolarV(:,:,3,k));
     hold on;


    
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Angular Velocity Comp Phidot Original'))
    grid on
    colorbar()
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    
    
    pause(1/fps)
end


%% Bar Graph Angular Velocity Filtered
figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Filtered Angular Velocity Model Bar Graph');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    b = bar3(polarV(:,:,2,k));
    hold on;
    
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Angular Velocity Comp Thetadot Filtered'))
    colorbar()
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    daspect([1 1 1])
    grid on
    %---------------------------------------------------------------------------------
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    b = bar3(polarV(:,:,3,k));
     hold on;


    
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Angular Velocity Comp Phidot Filtered'))
    grid on
    colorbar()
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    
    
    pause(1/fps)
end

%% Bar Graph Angular Acceleration Original
figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Original Angular Acceleration Model Bar Graph');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    b = bar3(originalpolarA(:,:,2,k));
    hold on;
    
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Angular Acceleration Comp Thetadot Original'))
    colorbar()
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    daspect([1 1 1])
    grid on
    %---------------------------------------------------------------------------------
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
   b = bar3(originalpolarA(:,:,3,k));
     hold on;


    
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Angular Acceleration Comp Phidot Original'))
    grid on
    colorbar()
    
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    
    pause(1/fps)
end
%% Bar Graph Angular Acceleration Filtered
figure('Units','normalized','Position',[0 0 1 1],'Name','Sail Filtered Angular Acceleration Model Bar Graph');  
hold on;
clf
kplot = 0;

startframe = 10;
endframe = 10;
numofrows = endframe - startframe + 1;
for k= startframe:endframe       %% Input frame interval that you want to see as limits for the for loop
    kplot = kplot + 1;
    
    
    subplot(numofrows,2,kplot)
    hold on;
%     plot3(Xpos(:,:,k),-Ypos(:,:,k),Zpos_flat(:,:,k),'.b');
    b = bar3(kalmanpolarA(:,:,2,k));
    hold on;
    
%     xlim([-1500 1750]);
%     ylim([-1500 1750]);
    %zlim([2500 3500]);
    title(strcat('Frame ',num2str(k),': Angular Acceleration Comp Thetadot Filtered'))
    colorbar()
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    daspect([1 1 1])
    grid on
    %---------------------------------------------------------------------------------
    kplot = kplot + 1;
    subplot(numofrows,2,kplot)
    hold on;
    b = bar3(kalmanpolarA(:,:,3,k));
     hold on;


    
    daspect([1 1 1])
    colormap(cmap(end:-1:1,:))
    title(strcat('Frame ', num2str(k),': Angular Acceleration Comp Phidot Filtered'))
    grid on
    colorbar()
    
    for i = 1:length(b)
    zdata = b(i).ZData;
    b(i).CData = zdata;
    b(i).FaceColor = 'interp';
    end
    
    pause(1/fps)
end