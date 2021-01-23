function [result1,result2,result3,result4]=fun_recon(sliceD,sliceB,Hdr,Param,reconType)
% Single slice reconstruction (CS+PI+PF) using different methods.

fsz=size(sliceD);

[lpB,hpB,lpD,hpD,hpS,lpS]=deal(zeros(size(sliceB)));

%%
switch reconType
    
    case 'KSPIC'         % K-space subtraction with intensity correction and phase correction (KSPIC)
        
        sliceS=sliceB-Hdr.IC*sliceD;    % Weighted subtraction
        
        for i=1:fsz(3)  % Acquire low-pass filtered centre and high-pass filtered data
            [lpB(:,:,i), hpB(:,:,i)] = fun_homo_filter(sliceB(:,:,i),Hdr,11);
            [lpD(:,:,i), hpD(:,:,i)] = fun_homo_filter(sliceD(:,:,i),Hdr,11);
            [lpS(:,:,i), hpS(:,:,i)] = fun_homo_filter(sliceS(:,:,i),Hdr,11);
        end
        
        IS_caliB = CScore_KSPIC(sliceB,hpS,Param,Param.sub);    % KSPIC reconstruction, using bright-blood data for SPIRiT calibration.
        IS_caliB_phaseD=fun_homo_phase(IS_caliB,lpD);           % Phase correction, using dark-blood data as phase reference.
        
        IS_caliB_phaseD_combB=fun_channel_combine('adpt-lowres',IS_caliB_phaseD,lpB);  % Adaptive channel combination, using bright-blood data as reference.
        result1=IS_caliB_phaseD_combB;
        
               
    case 'KS'         % Conventional K-space subtraction without intensity correction or phase correction (KS)
        
        sliceS=sliceB-sliceD;
        
        for i=1:fsz(3) 
            [lpS(:,:,i), hpS(:,:,i)] = fun_homo_filter(sliceS(:,:,i),Hdr,11);
        end
        
        IS = CScore_normal(hpS,Param,Param.sub);    % PI+CS reconstruction of subtracted data
        
        IS_homo=fun_homo_phase(IS,lpS);     % Partial Fourier homodyne reconstruction

        result1=fun_channel_combine('adpt-lowres',IS_homo,lpS); 
        
        
    case 'MS'       % Conventional magnitude subtraction with and without post intensity correction.  
        
        for i=1:fsz(3)     
            [lpD(:,:,i), hpD(:,:,i)] = fun_homo_filter(sliceD(:,:,i),Hdr,11);
            [lpB(:,:,i), hpB(:,:,i)] = fun_homo_filter(sliceB(:,:,i),Hdr,11);
        end
        
        IB = CScore_normal(hpB,Param,Param.mag);    % PI+CS reconstruction of bright-blood data
        ID = CScore_normal(hpD,Param,Param.mag);    % PI+CS reconstruction of dark-blood data
        IB_homo=fun_homo_phase(IB,lpB);             
        ID_homo=fun_homo_phase(ID,lpD);
        
        IB_homo_comb=fun_channel_combine('adpt-lowres',IB_homo,lpB);
        ID_homo_comb=fun_channel_combine('adpt-lowres',ID_homo,lpD);
        
        IS_homo_comb_IC=IB_homo_comb-Hdr.IC*ID_homo_comb;
        IS_homo_comb=IB_homo_comb-ID_homo_comb;
        
        result1=IS_homo_comb_IC; % MS-IC
        result2=IS_homo_comb;    % MS
        result3=IB_homo_comb;    % Bright-blood images
        result4=ID_homo_comb;    % Dark-bloodimages        
        
    otherwise
        error('Unsupport reconType');
end


end

