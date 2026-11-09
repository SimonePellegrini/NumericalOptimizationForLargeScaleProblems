function f = bandedTrigonometric_func(x)
   
    n = length(x);
    i_vec = (1:n)'; 
  
    sin_prev = [0; sin(x(1:n-1))];
    sin_next = [sin(x(2:n)); 0];
    
    terms = i_vec.*((1-cos(x))+sin_prev-sin_next);
    
    f = sum(terms);
end