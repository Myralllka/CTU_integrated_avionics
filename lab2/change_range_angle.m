function out = change_range_angle(angle,angle_last,angle_int,mode)
switch (mode)
    case 1 % range (0-2*pi) 
        angle = mod(angle,2*pi);
        if angle<0
            angle = angle + 2*pi;
        elseif angle >= 2*pi
            angle = angle - 2*pi;
        end
        out = angle;
     case 2 % without range in (rad) 
         delta = angle-angle_last;
         if (abs(delta) > pi)
             sign_delta = sign(delta);
             delta = (-1)*sign_delta*2*pi + angle - angle_last;
         end
         out = angle_int + delta;
end

end


       