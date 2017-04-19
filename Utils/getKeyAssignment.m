function b = getKeyAssignment(keyboard)KbName('UnifyKeyNames')% Function KbName converts between keycodes and key namesif keyboard==0 % Kits computer so windows keys    % keys for left side task    left1key = 'q';                            	% key for sure-absent response    left2key = 'w';                             % key for likely-absent response    left3key = 'e';                             % key for likely-present response    left4key = 'r';                           % key for sure-present response  "Clear" on some keyboards    b.respL = KbName({left1key,left2key,left3key,left4key});    % keys for right side task    right1key = 'u';                        	% key for sure-absent response    right2key = 'i';                             	% key for likely-absent response    right3key = 'o';                            	% key for likely-present response    right4key = 'p';                             	% key for sure-present response    b.respR = KbName({right1key,right2key,right3key,right4key});        b.abort = KbName('escape');else %lab computer    % keys for left side task    left1key = '1';                            	% key for sure-absent response    left2key = '4';                             % key for likely-absent response    left3key = '7';                             % key for likely-present response    left4key = 'tab';                           % key for sure-present response  "Clear" on some keyboards    b.respL = KbName({left1key,left2key,left3key,left4key});    % keys for right side task    right1key = 'enter';                        % key for sure-absent response    right2key = '+';                            % key for likely-absent response    right3key = '-';                            % key for likely-present response    right4key = 'delete';                       % key for sure-present response    b.respR = KbName({right1key,right2key,right3key,right4key});        b.abort = KbName('escape');end      otherKeys = KbName({'c','v','d','a','r','escape','space','return',...    'DownArrow','UpArrow'});  %for calibration, validation, drift correction% FlushEvents('keyDown');%Restrict keyboard so nothing works except necessary keysRestrictKeysForKbCheck([b.respL b.respR b.abort otherKeys]);