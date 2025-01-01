function s1(c)
N = sqrt(length(c));
hh = pcolor(reshape(c,N,N));

            set(hh, 'EdgeColor','none')
            shading interp
            hhh = colorbar; %caxis([0, 1]);
            % axis([0 2*pi pi-pi/gamma pi+pi/gamma])
            axis equal
            axis off
end