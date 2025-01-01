function writegif( filename , currentframe , framerate, withgcf )
% Write current figure into a gif movie
% input: filename, currentframe, framerate, doyouwantgcf

if withgcf ==1; frame = getframe(gcf); else; frame = getframe; end
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if currentframe == 1; imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'delaytime',1/framerate);
else; imwrite(imind,cm,filename,'gif','WriteMode','append','delaytime',1/framerate); end
end
