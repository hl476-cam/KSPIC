function [IC,Hdr]=fun_IC(rawD, rawB, Hdr, Param, interval,isImshow)
% To calculate the weighting factor for subtraction by using intensity correction based on robust regression.

fprintf('Intensity correction started.\n');
timeStart=clock;

fsz=Hdr.fullSize;

nSelect=1;

resultB_quick=zeros(fsz(2),fsz(3),ceil(fsz(1)/interval));
resultD_quick=zeros(fsz(2),fsz(3),ceil(fsz(1)/interval));

%%

% fprintf('Reconstructing slice:\r');

for iSlice=1:interval:fsz(1)
    
%     fprintf('%d\t',iSlice);
    
    sliceB=squeeze(rawB(iSlice,:,:,:));
    sliceD=squeeze(rawD(iSlice,:,:,:));
        
    [lpB,hpB,lpD,hpD]=deal(zeros(size(sliceB)));
    Hdr.centreSize=[12,6];
    for i=1:fsz(4)      % Acquire low-pass and high-pass data from subtracted data
        [lpD(:,:,i), hpD(:,:,i)] = fun_homo_filter(sliceD(:,:,i),Hdr,11);
        [lpB(:,:,i), hpB(:,:,i)] = fun_homo_filter(sliceB(:,:,i),Hdr,11);
    end
    
    [ID,IB] = CScore_quick(hpD,hpB);
    
    IB_pr=fun_homo_phase(IB,lpB);
    ID_pr=fun_homo_phase(ID,lpD);
    IB_pr_rss=fun_channel_combine('Psum',IB_pr,[]);
    ID_pr_rss=fun_channel_combine('Psum',ID_pr,[]);
    
    resultB_quick(:,:,nSelect)=IB_pr_rss;
    resultD_quick(:,:,nSelect)=ID_pr_rss;
    
    nSelect=nSelect+1;
    
end

%% IC based on robust regression

resultArrayB=abs(resultB_quick(:));
resultArrayD=abs(resultD_quick(:));

% tic
IC=robustfit_IC(resultArrayD(:),resultArrayB(:),'MRA5',Param.IC.alpha,[],'off');   

timeEnd=clock;
timeTotal=etime(timeEnd,timeStart);
fprintf('\nTotal intensity correction time=%dmin%ds\n\r\r',floor(timeTotal/60),round(rem(timeTotal,60)));

%% Plot
if isImshow
    col=dscatter_col(resultArrayD,resultArrayB);
    SCALE=0.6;
    BAR=0.5;
    msize = 2;
    marker = 'o';
    % set(figure,'position',[50 50 600 800]);
    figure;
    h = scatter(resultArrayD,resultArrayB,msize,col.^SCALE,marker,'filled');
    M=max([resultArrayD(:);resultArrayB(:)]);
    xlim([0 M]);ylim([0 M]);
    % xticks(0:2500:M);yticks(0:2500:M);
    hold on;plot([0 M],[0 IC*M],'--r','LineWidth',2);
    % plot([0 M],[0 Welsch*M],'--g','LineWidth',2);
    plot([0 M],[0 M],'--','LineWidth',2,'Color',[0 176/255 80/255]);
    hold off;
    colormap(jet);
    caxis([0 BAR]);
    axis square;
    title('Scatter map with density','FontSize',26,'FontWeight','bold');
    xlabel('Intensity on Dark-blood Images','FontSize',20,'FontWeight','bold');ylabel('Intensity on Bright-blood Images','FontSize',16,'FontWeight','bold');
    set(gca,'FontSize',16);
    set(gca,'LineWidth',2);
    set(gcf,'outerposition',get(0,'screensize'));
    
end

end


