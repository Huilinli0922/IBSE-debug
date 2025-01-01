function theta = hard_triangle(s)
theta = pi/2*(s<1/6) - pi/6*(s<1/2).*(s>=1/6) - pi*5/6.*(s>=1/2).*(s<5/6)...
            - 3*pi/2.*(s>=5/6);
theta(end/2+1) = (theta(end/2)+theta(end/2+2))/2;
        
  