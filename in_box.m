function box_num = in_box(pos, boxes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    box_num = 0;
    for i=1:size(boxes,1)
        if(pos(1) > boxes(i,1) && pos(1) < boxes(i,2) && pos(2) > boxes(i,3) && pos(2) < boxes(i,4))
            box_num = i;
            return;
        end
    end
end

