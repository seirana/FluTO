function [Ain, bin, tr_fc_lst, fctable] = add_new_tradeoffs(fctable, rep, Ain, bin, tr_fc_lst, m, f2, v2)
if rep > 1
    for j = 1:rep-1
        if tr_fc_lst(2,j) == tr_fc_lst(2,rep)
            eq = zeros(size(tr_fc_lst,1),1);
            A = zeros(size(tr_fc_lst,1),1);
            B = zeros(size(tr_fc_lst,1),1);
            
            for i = 3:size(tr_fc_lst,1)
                if tr_fc_lst(i,j) > 0 && tr_fc_lst(i,rep) > 0
                    eq(i,1) = tr_fc_lst(i,j);
                end
            end
            eq = unique(eq,'rows');
            if size(eq,1) == tr_fc_lst(2,rep)
                for i = 3:size(tr_fc_lst,1)
                    if tr_fc_lst(i,j) > 0 && tr_fc_lst(i,rep) == 0
                        A(i,1) = i;
                    end
                    if tr_fc_lst(i,j) == 0 && tr_fc_lst(i,rep) > 0
                        B(i,1) = i;
                    end
                end
                
                A = unique(A,'rows');
                B = unique(B,'rows');
                
                for li = 2:size(A,1)
                    for ri = 2:size(B,1)
                        fctable(A(li,1)-2,B(ri,1)-2) = 5;
                        fctable(B(ri,1)-2,A(li,1)-2) = 5;
                        s_in = zeros(1,v2);
                        s_in(1,A(li,1)-2) = 1;
                        s_in(1,B(ri,1)-2) = 1;
                        new_bin = [zeros(1,m), zeros(1,f2+v2), zeros(1,v2), s_in];
                        Ain = [Ain; new_bin];
                        bin = [bin; 1];
                    end
                end
            end
        end
    end
end

end
