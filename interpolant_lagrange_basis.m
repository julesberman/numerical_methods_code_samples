function Pofx = evalInterpolantLagrangeBasis(xArr, t_nodes, f_nodes)
	Pofx = zeros(size(xArr));
	for j = 1:length(xArr)
        x = xArr(j);
	    for k = 1:length(t_nodes)
            LkofX = ones(size(xArr));
	        for i = 1:length(t_nodes)
	            if i ~= k

	                xi = t_nodes(i);
	                xk = t_nodes(k);
	                LkofX(j) = LkofX(j) * ((x - xi) / (xk - xi));
	            end 
	        end
	        yk = f_nodes(k);
	        Pofx(j) = Pofx(j) + (LkofX(j)*yk);
	    end
	end
end