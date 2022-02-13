filepath = "Earth_bottom";

load(filepath + "/Xraw.mat");
load(filepath + "/Yraw.mat");
load(filepath + "/Zraw.mat");
load(filepath + "/Ampraw.mat");

[pxX,pxY,nbrOfFrames] = size(Ampraw);

nbrOfFrames2 = floor(nbrOfFrames / 10) + min(mod(nbrOfFrames,10),1);

Xraw2 = zeros(pxX,pxY,nbrOfFrames2);
Yraw2 = zeros(pxX,pxY,nbrOfFrames2);
Zraw2 = zeros(pxX,pxY,nbrOfFrames2);
Ampraw2 = zeros(pxX,pxY,nbrOfFrames2);


%%% We get rid of the values where Zraw=-1 or 0 %%%
index = (Zraw==-1);
Xraw(index) = nan;          
Yraw(index) = nan;        
Zraw(index) = nan;       
Ampraw(index) = nan;     
index = (Zraw==0);
Xraw(index) = nan;         
Yraw(index) = nan;       
Zraw(index) = nan;          
Ampraw(index) = nan;       


for k=1:nbrOfFrames2
    
    Beg = (k-1)*10+1;
    End = min(k*10,nbrOfFrames);
    
    Xraw2(:,:,k) = mean_global_matrix(Xraw(:,:,Beg:End),numel(Beg:End));
    Yraw2(:,:,k) = mean_global_matrix(Yraw(:,:,Beg:End),numel(Beg:End));
    Zraw2(:,:,k) = mean_global_matrix(Zraw(:,:,Beg:End),numel(Beg:End));
    Ampraw2(:,:,k) = mean_global_matrix(Ampraw(:,:,Beg:End),numel(Beg:End));
    
end

Xraw_avg = Xraw2;
Yraw_avg = Yraw2;
Zraw_avg = Zraw2;
Ampraw_avg = Ampraw2;

%%% We get rid of the values where Zraw=nan %%%
index = (isnan(Zraw_avg));
Xraw_avg(index) = -1;          
Yraw_avg(index) = -1;        
Zraw_avg(index) = -1;       
Ampraw_avg(index) = -1;   


save(filepath + "\Xraw_avg.mat", 'Xraw_avg');
save(filepath + "\Yraw_avg.mat", 'Yraw_avg');
save(filepath + "\Zraw_avg.mat", 'Zraw_avg');
save(filepath + "\Ampraw_avg.mat", 'Ampraw_avg');