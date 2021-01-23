function [Ic_real,Ic] = fun_homo_phase(I,lpK)
% Homodyne phase correction by using the central region of k-space;

lpI=fftshift(flip(flip(ifft2(fftshift(fftshift(lpK,1),2)),1),2),1);   % Do 2D fft to obtain the low-pass image

Ic=I./exp(1i*angle(lpI));   % Phase correction
Ic_real=real(Ic);
% Ic(Ic<0)=0;

end