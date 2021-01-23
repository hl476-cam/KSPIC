function I_combined=fun_channel_combine(mode,I_channel,K_ref_low)

% Channel combination using adaptive combination or root-of-squares

switch mode
    case 'adpt-lowres'
        % Method1: adaptive combination using a low resolution dark-blood image
        I_ref=fftshift(flip(flip(ifft2(fftshift(fftshift(K_ref_low,1),2)),1),2),1);   % Do 2D fft to obtain the low-pass image
        I_ref_c=I_ref./exp(1i*angle(I_ref));   % Homodyne correction / Phase recovery
        I_ref_c=real(I_ref_c);
        I_combined=adapt_array_2d(I_ref_c,I_channel);
               
    case 'Psum'
        % Method2: root-of-squares with polarity preserved.
        Psum=sum(sign(I_channel).*abs(squeeze(I_channel)).^2,3);
        I_combined=sign(Psum).*sqrt(abs(Psum));
        
    case 'rss'
        % Method3: root-of-squares
        I_channel(I_channel<0)=0;
        I_combined=sqrt(sum(abs(squeeze(I_channel)).^2,3));
        
end

I_combined=flip(I_combined,1);

end

