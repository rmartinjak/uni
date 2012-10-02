source ab2_kmeans.m

function [L, Lvec, i] = em (Xmat, K)
	[d, N] = size(Xmat);

	% initialize Umat, Dmat, Hmat with k-means
	[_, _, _, _, Umat] = kmeans(Xmat, K, @RandChoiceProto);
	Dmat = mkDist(Xmat, Umat);
	Hmat = mkAssign(Dmat);
	sigmaSquare = sum(sum(Dmat' .* Hmat)) / (N*d);

	% P(j)
	Pmat = sum(Hmat) ./ N;

	% likelihood vector
	Lvec = zeros(1, ITER_MAX);
	Lvec(1) = inf;

	i = 1;
	while ((i == 1) || abs(Lvec(i) - Lvec(i-1)) > EPSILON) && ++i < ITER_MAX

		% p(x|j)
		pmat = (2 * pi * sigmaSquare) ** -(d/2) * exp( (-1 / (2*sigmaSquare)) * Dmat' );

		% assignment matrix h(i,j)
		Hmat = 	(repmat(Pmat, rows(pmat), 1) .* pmat) ./ ... % p(x|j)*P(j)	
				repmat(
					sum( repmat(Pmat, rows(pmat), 1) .* pmat, 2 ),
					1,
					columns(pmat)
				);

		% P(j)
		Pmat = sum(Hmat) ./ N;

		% µ(j)
		Umat = Xmat * Hmat ./ repmat(sum(Hmat), rows(Xmat), 1);

		% x_i - mu_j
		Dmat = mkDist(Xmat, Umat);

		sigmaSquare = sum(sum(Hmat .* Dmat')) / (N*d);

		Lvec(i) = sum(
					log(
						sum(repmat(Pmat, rows(pmat), 1) .* pmat, 2)	% p(x_i)
					)
				);
		L = Lvec(i);
	end

end



function sheet2_ex1a(Xmat)
	[_, Lvec, Iter] = em(Xmat, 4);
	plot([2:Iter], Lvec(2:Iter))
end

function sheet2_ex1b(Xmat, Kmax)
	Lvec = zeros(1, Kmax);
	for K = 2:Kmax
		[L, _, _] = em(Xmat, K);
		Lvec(K) = L;
		L
	end
	plot([2 : Kmax], Lvec(2 : end));
end


function sheet2_ex2(Xmat, Kmax)
	[d, N] = size(Xmat);
	BICvec = zeros(1, Kmax);
	for K = 2:Kmax
		[L, _, _] = em(Xmat, K);
		M = K*d + (K-1) + 1;
		BICvec(K) = 2 * L - M * log(N);
	end
	plot([2 : Kmax], BICvec(2 : end));
end
	
%############%
% Exercise 1 %
%############%
% sheet2_ex_1a(GenXmat(0.25));

% load GexprData1.mat;
% sheet2_ex_1b(Xmat, 10);

%############%
% Exercise 2 %
%############%
% load GexprData1.mat;
% sheet2_ex_2(Xmat, 10);
