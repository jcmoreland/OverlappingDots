function drawFix(display)
%drawFix.m
%makes fixation cross

if ~isfield(display,'fixColor')
    display.fixColor = 255*ones(1,3);
end

if ~isfield(display,'fixOffset')
    display.fixOffset = [0 0];
end
center = display.resolution/2 + angle2pix(display,display.fixOffset);

%horizontal line
Screen('DrawLine', display.windowPtr, display.fixColor, center(1), center(2)-8, center(1), center(2)+8, 3)
%vertical line
Screen('DrawLine', display.windowPtr, display.fixColor, center(1)-8, center(2), center(1)+8, center(2), 3)
