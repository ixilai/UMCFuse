function rp = rotatePosition(w)
    p = (1:w)';
    rp = zeros(w, 722);
    for a = 0:0.5:180
        if a<=90
            b = a;
        end
        if a>90
            b = a-180;
        end
        px = floor( p*cosd(b)+0.5 ); py = floor( p*sind(b)+0.5 );
        rp(:, 2*a+1) = px; rp(:, 2*a+362) = py;
    end
end

