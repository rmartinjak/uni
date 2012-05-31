
# data vectors
Dmat = [ randn(2, 100).*0.5, \
	randn(2, 100).*0.5 + [0; 1], \
	randn(2, 100).*0.5 + [1; 0], \
	randn(2, 100).*0.5 + [1; 1] ]

K = 4

Pmat = rand(2, K)


# distance matrix
for i = 1:size(v)(:,1)
	for k = 1:size(p)(:,1)
		D(i,k) = sumsq(v(i,:) - p(k,:))
	end
end

