clear all; 
close all;
%input files 
fIn1 = "Q:\Users\yc645\20250320_1_turbulent_100g_clockwise_subregion.dfm";
fIn2 = "Q:\Users\yc645\20250320_background_subregion.dfi";
fIn3 = fIn2;
fIn4 = "Q:\Users\yc645\20250320_1_uniform_concentration_subregion.dfi";
fIn5 = fIn4;
x_map_file = "Q:\Users\yc645\LIF-processing-main\LIF-processing-main\digiCode\maybeGood_x_map.txt";
y_map_file = "Q:\Users\yc645\LIF-processing-main\LIF-processing-main\digiCode\maybeGood_y_map.txt";

%read input files
backBrightDfi = dfi2mat(fIn2);
backDarkDfi = dfi2mat(fIn3);
DyeBrightDfi = dfi2mat(fIn4);
DyeDarkDfi = dfi2mat(fIn5);
blackBright = backBrightDfi.image;
blackDark = backDarkDfi.image;
DyeBright = DyeBrightDfi.image;
DyeDark = DyeDarkDfi.image;

%Concentration of unifrom dye images 
C0 = 0.5;

%Domain Details
yMin = 0;
yMax = 2047;
xMin = 0; 
xMax = 1731;
w = xMax-xMin;
h = yMax-yMin;

%subtranct black from background 
backBright = DyeBright - blackBright;
backDark = DyeDark - blackDark;

tic
%initialise buffer
buffLen = 50;
startFrame = 100;
endFrame = 101;
nBuff = ceil((endFrame - startFrame)/buffLen);
p1 = df_dfm_info(fIn1);
f1 = fopen(fIn1,'r');
%read in buffLen frames from input movie  
for n = 1:nBuff
    pos = startFrame + (n-1)*buffLen;
    if (((endFrame-startFrame) - n*buffLen)/buffLen > 1) 
        frms1=pos:(pos+buffLen-1);
    else 
        frms1=pos:endFrame;
    end 
    exp = df_dfm_read(f1, frms1, p1);
    Len = length(frms1);
%initialise parallel for loop to process frames in buffer
    for ii = 1:Len
        i = (ii-1) + pos;
        im = dfiRegion(double(exp(:,:,1,ii))/255,xMin,xMax,yMin,yMax);
        if (mod(i,2)==0)  
            rhoI = double(im) - blackBright;
            back = backBright;
            rhoI(rhoI<0) = 0;
            rhoI = filtWindow(rhoI,15,1,2);
        else 
            rhoI = double(im) - blackDark;
            back = backDark;
            rhoI(rhoI<0) = 0;
            rhoI = filtWindow(rhoI,15,1,2);
        end
        rhoIRaw = rhoI;
        %corr(raw plif, background, x map to laser sheet space, y map (equiv), width, height, have streak, want synthetic laser sheet)
        %corr = correct_Plif_MI(rhoI,back,x_map_file,y_map_file,w,h,C0,true,true,200,5,1,8,4);
        [corr,mapped,streakRM,illFix] = correct_Plif_MI_new(rhoI,back,x_map_file,y_map_file,w,h,C0,true,true,200,20,1,10,4);
        %imwrite(corr,"..\outputs\strat[20250317-3]_"+i+".tif");
    end
end
toc
fclose('all');

X = linspace(0,w-1,w);
Y = linspace(h-1,0,h);
figure;
imagesc(X,Y,corr);
colormap(flipud(bone));
colorbar;
caxis([0 1]);
pbaspect([1 h/w 1]);
axis xy; 
axis off;
ax=gca;
%exportgraphics(ax,"corr.jpg",'BackgroundColor','none')

X = linspace(0,w-1,w);
Y = linspace(h-1,0,h);
figure;
imagesc(X,Y,rhoIRaw);
colormap(flipud(jet));
%colorbar;
caxis([0 1]);
pbaspect([1 h/w 1]);
axis xy; 
axis off;
ax=gca;
%exportgraphics(ax,"input.jpg",'BackgroundColor','none')

X = linspace(0,w-1,w);
Y = linspace(h-1,0,h);
figure;
imagesc(X,Y,flipud(mapped));
colormap(flipud(jet));
%colorbar;
%caxis([0 1]);
pbaspect([1 h/w 1]);
axis xy; 
axis off;
ax=gca;
%exportgraphics(ax,"inputMapped.jpg",'BackgroundColor','none')

X = linspace(0,w-1,w);
Y = linspace(h-1,0,h);
figure;
imagesc(X,Y,streakRM);
colormap(flipud(jet));
%colorbar;
caxis([0 1]);
pbaspect([1 h/w 1]);
axis xy; 
axis off;
ax=gca;
%exportgraphics(ax,"StreakRm.jpg",'BackgroundColor','none')

X = linspace(0,w-1,w);
Y = linspace(h-1,0,h);
figure;
imagesc(X,Y,illFix);
colormap(flipud(jet));
colorbar;
caxis([0 1]);
pbaspect([1 h/w 1]);
axis xy; 
axis off;
ax=gca;
%exportgraphics(ax,"corrMapped.jpg",'BackgroundColor','none')

rhoProf = mean(corr(:,round(2*w/5):round(3*w/5)),2);
z = linspace(0,h-1,h);
figure;
h1 = axes;
plot(fliplr(rhoProf),z);
set(h1,'Ydir', 'reverse')



