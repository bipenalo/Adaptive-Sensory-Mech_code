% function to estimate the neuron location at which activity is equal to th

function sp = EstimateCell(inp, t, th, space, ch) 

    %space_value1 = linspace(1,1800,1800);
    S = inp(t, space(1):space(2)); % get the neuron location at 70 units of activity
    x = find(S > th);
    
    % get the interpolated value of x given y
    %sp = interp1(space8_1, space_value1, th);
    if isempty(x)
        sp = 0;
    elseif ch == 1
        sp = x(1);
    else
        sp = 1300 + x(end); % add the initial location at 1300 to the final element which was bigger than th!
    end
end
