% the evaluate the hermite interpolation directly
function Pofx = evalInterpolantHermite(xArr, t_nodes, f_nodes, df_nodes)
	Pofx = zeros(size(xArr));
	for j = 1:length(xArr)
        x = xArr(j);
	    for k = 1:length(t_nodes)
            
            LkofX = ones(size(xArr));
            HkofX = ones(size(xArr));
            KkofX = ones(size(xArr));
            xk = t_nodes(k);
            dLkofxk = 0;
	        for i = 1:length(t_nodes)
	            if i ~= k
	                xi = t_nodes(i);
	                LkofX(j) = LkofX(j) * ((x - xi) / (xk - xi));
                    dLkofxk = dLkofxk + 1/(xk - xi);
                end 
            end
     
            HkofX(j) = LkofX(j)^2*(1-2*dLkofxk*(x - xk));     
            KkofX(j) = LkofX(j)^2*(x - xk);
                    
	        yk = f_nodes(k);
            zk = df_nodes(k);
	        Pofx(j) = Pofx(j) + (HkofX(j)*yk)+ (KkofX(j)*zk);
	    end
	end
end