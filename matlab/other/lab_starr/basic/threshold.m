function y = threshold(x,t,d)
% function returns matrix of indices at which the 1-dimensional 
% matrix y crosses threshold t.  The direction is from positive 
% to negative values if d = '+-', and from negative to positive,
% if d = '-+'.  If d is '+-+', the function will return a 
% 1-dimensional matrix with start- and end indices for +-+ 
% threshold crossings.  The reverse is true for '-+-'.  The 
% function will return a value of -1 if no threshold crossings 
% are identified.
%
% Example:
%		a = threshold(M,5,'-+')
%
% This will return a matrix of numbers representing the indices 
% at which a value of 5 was crossed from "below"
% 
y = -1;							% sets default return value
if length(d) == 3
   if d == '+-+'
      for i = 2:length(x)     % for loop evaluates neighboring 
      	if x(i) < t 			% values for threshold criterion + -> -
      		if x(i-1) >= t
	            y = [y i];
   	      end
    		end
   	end
   	if length(y) > 2
     	 y1 = y(2);
      	for i = 2:length(x)     % for loop evaluates neighboring 
   	   	if x(i) >= t 			% values for threshold criterion - -> +
      			if x(i-1) < t
            	y = [y i-1];
         		end
  		  		end
   		end
   		y = sort(y);
   		if y(1) ~= y1
      		 y = y(2:length(y));
   		end
   		if rem(length(y),2)~=1
      		 y = y(1:length(y)-1);
       	end
   	end
	elseif d == '-+-'
      for i = 2:length(x)     % for loop evaluates neighboring
      	if x(i) > t 			% values for threshold criterion - -> +
      		if x(i-1) <= t
               y = [y i];
         	end
    		end
    	end
    	if length(y) > 1
    		y1 = y(2);
    		for i = 2:length(x)     % for loop evaluates neighboring 
      		if x(i) <= t 			% values for threshold criterion + -> +
      			if x(i-1) > t
            		y = [y i-1];
         		end
    			end
    		end
    		y = sort(y);
    		if y(2) ~= y1
      		y = [-1 y(3:length(y))];
    		end
    		if rem(length(y),2)~=1
      		y = y(1:length(y)-1);
      	end
      end 
   end
else if length(d) == 2  
	if d == '-+'					% if loop takes care of direction 
  		x = -x;						% of threshold crossing
   	t = -t;
	end
	for i = 2:length(x)			% for loop evaluates neighboring 
   	if x(i) < t 				% values for threshold criterion
      	if x(i-1) >= t
         	y = [y i];
      	end
   	end
   end
  end
end
if length(y) > 1				% if loop removes the default return
   y = y(2:length(y));		% value if threshold crossings are found
end
