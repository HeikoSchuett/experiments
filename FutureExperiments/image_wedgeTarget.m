function imageBitmap = image_wedgeTarget(imSize,nWedges,propFilled,pxTarget,pxSurround1,pxSurround2,pxBackground,targetC,surroundC,bgC,wedgeC)
%function imageBitmap = image_wedgeTarget(imSize,pxTarget,pxSurround,pxBackground,targetC,surroundC,wedgeC)
% this function creates a wedgeTarget stimulus with a certain imSize in
% pixels.
% the arguments pxXXXX give the radii of the circles in the stimulus.
% targetC and surroundC give the gray values for the target & surround
% wedgeC should be a 4 element vector with values for the 4 wedges

assert(propFilled<=1,'proportion filled must be <= 1');
assert(propFilled>=0,'proportion filled must be >= 0');

antialiasing = 8;

imSize = antialiasing*imSize;

imageBitmap = bgC*ones(imSize);
x = 1:imSize(2);
y = 1:imSize(1);
x = x-floor(imSize(2)/2);
y = y-floor(imSize(1)/2);

x = x./antialiasing;
y = y./antialiasing;

[xx,yy]=meshgrid(x,y');

r = sqrt(xx.^2+yy.^2);
phi = atan2(yy,xx);
phi2 = phi;
phi2(phi2<0) = phi2(phi2<0)+2*pi;

target = r<pxTarget;
surround = (r>=pxTarget)&(r<pxSurround1);
background = (r>=pxSurround2).*(r<pxBackground)+(r>=pxSurround1).*(r<pxSurround2).*(0.5-0.5*cos(pi*(r-pxSurround1)./(pxSurround2-pxSurround1)));

wedges = zeros(size(imageBitmap));
wAngles = linspace(-pi,pi,nWedges+1);
anglePerWedge = 4*pi/nWedges*propFilled;
for iw = 1:nWedges
    if (wAngles(iw)+anglePerWedge)<=pi
        wedges = wedges+(background.*(phi>=wAngles(iw)).*(phi<(wAngles(iw)+anglePerWedge))...
            .* (0.5-0.5*cos(2.*pi./anglePerWedge.*(phi-wAngles(iw)))));
    else
        wedges = wedges+(background.*(phi2>=wAngles(iw)).*(phi2<(wAngles(iw)+anglePerWedge))...
            .* (0.5-0.5*cos(2.*pi./anglePerWedge.*(phi2-wAngles(iw)))));
    end
    
end

imageBitmap(target) = targetC;
imageBitmap(surround) = surroundC;
imageBitmap = imageBitmap+wedges.*wedgeC;

imageBitmap = imresize(imageBitmap,1/antialiasing,'box');