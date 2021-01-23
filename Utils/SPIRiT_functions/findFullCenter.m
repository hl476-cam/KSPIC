function Kcalib = findFullCenter(SampleMask)
%init
xsize=2;
ysize=2;

SMsz=size(SampleMask);

Kcalib=mycrop(SampleMask,[xsize,ysize],'c');

area=prod([xsize,ysize]);
isfullysampled=(sum(Kcalib(:))==area);

%grow a square
while(isfullysampled && (xsize<SMsz(1)) && (ysize<SMsz(2)) )
    xsize=xsize+1;
    ysize=ysize+1;
    
    Kcalib=mycrop(SampleMask,[xsize,ysize],'c');
    
    area=prod([xsize,ysize]);
%     disp([num2str(xsize) '^2 sum(Kcalib(:)) ',num2str(sum(Kcalib(:))),' == ', num2str(area),' area']);
    isfullysampled=(sum(Kcalib(:))==area);
end

%go back one step
xsize=xsize-1;
ysize=ysize-1;
Kcalib=mycrop(SampleMask,[xsize,ysize],'c');
    
isfullysampled=1;
%try growing a rectangle
if(SMsz(1)>=SMsz(2))
    while(isfullysampled && xsize<SMsz(1))
        xsize=xsize+1;
%         disp([num2str(xsize) 'x' num2str(ysize)]);
        Kcalib=mycrop(SampleMask,[xsize,ysize],'c');
        
        area=prod([xsize,ysize]);
        isfullysampled=(sum(Kcalib(:))==area);
    end
    
    %go back one step
    xsize=xsize-1;
    Kcalib=mycrop(SampleMask,[xsize,ysize],'c');
else
    while(isfullysampled && ysize<SMsz(2))
        ysize=ysize+1;
        
        Kcalib=mycrop(SampleMask,[xsize,ysize],'c');
        
        area=prod([xsize,ysize]);
        isfullysampled=(sum(Kcalib(:))==area);
    end
    
    %go back one step
    ysize=ysize-1;
    Kcalib=mycrop(SampleMask,[xsize,ysize],'c');
end

return;

    
    