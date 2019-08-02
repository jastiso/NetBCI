function [] = check_error(err,tol)
%compare error to tolerance, print warning

if err > tol
  warning ('Your error is too large. Check that you are reaching your target state. If not, try adding more entries to B, fewer to S, or using smaller matrices')
end

end

