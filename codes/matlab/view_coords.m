function view_coords(x)
    boundary = 5;
    c = mean(x, 2);
    for i = 1:length(c)
        x(i,:) = x(i,:) - c(i);
    end
    for i =	1:length(x)
        plot3(x(1,i), x(2,i), x(3,i), '*b');
        if i == 1
            hold on
            box  on;
            grid on;
            xlim([-1,1] * boundary)
            ylim([-1,1] * boundary)
            zlim([-1,1] * boundary)
        end
    end
    hold off;
end