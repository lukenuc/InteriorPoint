% function [urow, ucol, uval, c, exitflag] = GaussElimCOO(row, col, val, b, n)
% 
% exitflag = 0;  % Success flag
% urow = row;
% ucol = col;
% uval = val;
% c = b;
% 
% 
% for j = 1:n
% 
%     pivot_idx = find(urow == j & ucol == j, 1);
% 
%     if isempty(pivot_idx)
%         exitflag = -1;
%         return
%     end
% 
%     ajj = uval(pivot_idx);  % Pivot element
% 
%     % Find elements below the pivot in column j
%     row_idx = find(ucol == j & urow > j);
% 
%     for idx = row_idx'
%         i = urow(idx);
%         aij = uval(idx);
%         m = aij / ajj;
% 
%         % Eliminate element
%         uval(idx) = 0;
% 
%         % Update RHS vector
%         c(i) = c(i) - m * c(j);
% 
%         % Process row i
%         temp = find(urow == j);
%         for cnt = 1:length(temp)  % Elements in pivot row
%             k = temp(cnt);
%             t = ucol(k);
%             if t == j
%                 continue  % Skip diagonal
%             end
% 
%             fetched_rows = urow == i; 
%             fetched_cols = ucol(fetched_rows); 
%             existing_idx = find(fetched_cols == t, 1); 
% 
%             if isempty(existing_idx)
%                 % If entry does not exist, insert new nonzero
%                 urow = [urow; i];
%                 ucol = [ucol; t];
%                 uval = [uval; -m * uval(k)];
%             else
%                 % Otherwise, update existing value
%                 uval(existing_idx) = uval(existing_idx) - m * uval(k);
%             end
%         end
%     end
% 
%     % % Remove zeroed elements (optional optimization)
%     % nonzero_mask = uval ~= 0;
%     % urow = urow(nonzero_mask);
%     % ucol = ucol(nonzero_mask);
%     % uval = uval(nonzero_mask);
% end
% 
% end


% function [urow, ucol, uval, c, exitflag] = GaussElimCOO(row, col, val, b, n)
% 
% exitflag = 0;
% k_ptr = 1;
% row_ptr = 1;
% col_ptr = 1;
% val_ptr = 1;
% j = 1;
% urow = row; ucol = col; uval = val; c = b; 
% 
% while j < n
% 
%     col_ptr = find(urow == j & ucol == j, 1);
%     if (~isempty(col_ptr))
%         ajj = uval(col_ptr);
%         col_ptr = col_ptr + 1;
%         while ucol(col_ptr) == j
%             row_ptr = col_ptr;
%             i = urow(row_ptr);
%             aij = uval(row_ptr);
%             m = aij/ajj;
%             uval(row_ptr) = []; urow(row_ptr) = []; ucol(row_ptr) = [];
%             c(i) = c(i) - m*c(j); 
%             k_ptr = row_ptr + 1;
%             while k_ptr <= length(urow)
%                 if urow(k_ptr) == i
%                     t_ptr = k_ptr; 
%                     t = ucol(k_ptr); 
%                     while ucol(t_ptr) == t
%                         if urow(t_ptr) == j 
%                             uval(k_ptr) = uval(k_ptr) - m*uval(t_ptr); 
%                             break
%                         end
%                         t_ptr = t_ptr - 1; 
%                     end
%                 end
%                 k_ptr = k_ptr + 1;
%             end
%         end
%         j = j+1; 
%     else
%         exitflag = -1;
%         break
%     end
% 
% end
% 
% end

function [urow, ucol, uval, c, exitflag] = GaussElimCOO(row, col, val, b, n)

exitflag = 0;
k_ptr = 1;
row_ptr = 1;
col_ptr = 1;
val_ptr = 1;
j = 1;
urow = row; ucol = col; uval = val; c = b; 

while j < n

    col_ptr = find(urow == j & ucol == j, 1);
    if (~isempty(col_ptr))
        ajj = uval(col_ptr);
        col_ptr = col_ptr + 1;
        while ucol(col_ptr) == j
            row_ptr = col_ptr;
            i = urow(row_ptr);
            aij = uval(row_ptr);
            m = aij / ajj;
            
            % Zero out the element at (i, j) and update right-hand side
            uval(row_ptr) = []; urow(row_ptr) = []; ucol(row_ptr) = [];
            c(i) = c(i) - m * c(j); 

            k_ptr = row_ptr + 1;
            while k_ptr <= length(urow)
                if urow(k_ptr) == i
                    t_ptr = k_ptr; 
                    t = ucol(k_ptr); 
                    while ucol(t_ptr) == t
                        if urow(t_ptr) == j 
                            % Update the value at (i, t) based on the multiplier
                            uval(k_ptr) = uval(k_ptr) - m * uval(t_ptr); 
                            
                            % Add a check for when an element becomes small or zero
                            % Apply a tolerance for floating point precision
                            if abs(uval(k_ptr)) < 1e-10
                                uval(k_ptr) = 0;
                            end
                            break;
                        end
                        t_ptr = t_ptr - 1; 
                    end
                end
                k_ptr = k_ptr + 1;
            end
        end
        j = j + 1; 
    else
        exitflag = -1;
        break
    end

end

end
