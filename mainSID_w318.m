clear all; 
close all;

% 
pool = gcp('nocreate'); % Returns the current pool object or empty if none exists
if isempty(pool)
    disp('No active parallel pool.');
    numWorkers = 10;
    parpool('local', numWorkers);
else
    disp(['Parallel pool with ', num2str(pool.NumWorkers), ' workers is active and will be closed']);
    delete(gcp('nocreate'));
    numWorkers = 10;
    parpool('local', numWorkers);
end

%input files 
fIn1 = "V:\202311\w318\LIF-Processing\inputs\CamC_dimmer.dfm";
fIn2 = "V:\202311\w318\LIF-Processing\inputs\CamC_black.dfi";
fIn3 = fIn2;
fIn4 = "V:\202311\w318\LIF-Processing\inputs\CamC_bright.dfi";
fIn5 = fIn4;
x_map_file = "V:\202311\w318\LIF-Processing\inputs\w318_ribbon_test2_x_map.txt";
y_map_file = "V:\202311\w318\LIF-Processing\inputs\w318_ribbon_test2_y_map.txt";

x_map_lab = "V:\202311\w318\LIF-Processing\inputs\w318_target_x_map.txt";
y_map_lab = "V:\202311\w318\LIF-Processing\inputs\w318_target_y_map.txt";

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
yMax = 1024;
xMin = 0; 
xMax = 3320;
w = xMax-xMin;
h = yMax-yMin;

%subtranct black from background 
backBright = DyeBright - blackBright;
backDark = DyeDark - blackDark;

tic
%initialise buffer
buffLen = 50;
startFrame = 49000;
endFrame =50200;
nBuff = ceil((endFrame - startFrame + 1) / buffLen);
corr_stack = zeros([h, w, endFrame-startFrame +1]);
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
    frameIndices = frms1;
    parfor ii = 1:Len
        frameIndex = frameIndices(ii);
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
        % corr = correct_Plif_MI(rhoI,back,x_map_file,y_map_file,w,h,C0,true,true,200,20,1,10,4);
        %[corr,mapped,streakRM,illFix] = correct_Plif_MI_new(rhoI,back,x_map_file,y_map_file,w,h,C0,true,true,200,20,1,10,4);
        % imwrite(corr,"..\outputs\strat[20250318_1_t_50]_"+i+".tif");
        %stackInd = i - startFrame + 1;
        %corr_stack(:,:,ii-startFrame+pos) = corr;
        corr_stripe = correct_Plif_MI_20250906(rhoI,back,x_map_file,y_map_file,w,h,C0,true,true,200,15,1,15,4);
        window_size = 15;     
        filtType    = 1;     % Stefan: normally use type 1
        corr_filt = filtWindow(corr_stripe, window_size, filtType);% written by star
        corr = mapTo(corr_filt, x_map_lab, y_map_lab, h, w, false);      
        imwrite(corr, sprintf('..\\outputs\\w318_test5_%05d.tif', i));
        %filename = sprintf('V:\202311\w314\LIF-processing-main\LIF-processing-main\outputs\strat[20251201_1_trail]_%05d.tif', frameIndex);
        %imwrite(corr, filename);
    end
end
toc
fclose('all');

%save('corr_0306c50.mat', 'corr_stack', '-v7.3');


% X = linspace(0,w-1,w);
% Y = linspace(h-1,0,h);
% figure;
% imagesc(X,Y,corr);
% colormap(flipud(bone));
% colorbar;
% %caxis([0 1]);
% pbaspect([1 h/w 1]);
% axis xy; 
% axis off;
% ax=gca;
% exportgraphics(ax,"corr.jpg",'BackgroundColor','none')
% 
% X = linspace(0,w-1,w);
% Y = linspace(h-1,0,h);
% figure;
% imagesc(X,Y,rhoIRaw);
% colormap(flipud(jet));
% %colorbar;
% %caxis([0 1]);
% pbaspect([1 h/w 1]);
% axis xy; 
% axis off;
% ax=gca;
% exportgraphics(ax,"input.jpg",'BackgroundColor','none')
% 
% % X = linspace(0,w-1,w);
% % Y = linspace(h-1,0,h);
% % figure;
% % imagesc(X,Y,flipud(mapped));
% % colormap(flipud(jet));
% % %colorbar;
% % %caxis([0 1]);
% % pbaspect([1 h/w 1]);
% % axis xy; 
% % axis off;
% % ax=gca;
% % %exportgraphics(ax,"inputMapped.jpg",'BackgroundColor','none')
% 
% X = linspace(0,w-1,w);
% Y = linspace(h-1,0,h);
% figure;
% imagesc(X,Y,streakRM);
% colormap(flipud(jet));
% %colorbar;
% caxis([0 1]);
% pbaspect([1 h/w 1]);
% axis xy; 
% axis off;
% ax=gca;
% %exportgraphics(ax,"StreakRm.jpg",'BackgroundColor','none')
% 
% X = linspace(0,w-1,w);
% Y = linspace(h-1,0,h);
% figure;
% imagesc(X,Y,illFix);
% colormap(flipud(jet));
% colorbar;
% caxis([0 1]);
% pbaspect([1 h/w 1]);
% axis xy; 
% axis off;
% ax=gca;
% %exportgraphics(ax,"corrMapped.jpg",'BackgroundColor','none')
% 
% rhoProf = mean(corr(:,round(2*w/5):round(3*w/5)),2);
% z = linspace(0,h-1,h);
% figure;
% h1 = axes;
% plot(fliplr(rhoProf),z);
% set(h1,'Ydir', 'reverse')
% 
% delete(gcp('nocreate'));

