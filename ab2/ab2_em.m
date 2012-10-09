function ab2_em
%############%
% Exercise 1 %
%############%
% sheet2_ex1a(GenXmat(0.25));

% load GexprData1.mat;
% sheet2_ex1b(Xmat, 10);

%############%
% Exercise 2 %
%############%
% load GexprData1.mat;
% sheet2_ex2(Xmat, 10);
end

function Pmat = RandChoiceProto(Xmat, K)
    Perm = randperm(columns(Xmat));
	Pmat = Xmat(:,Perm(1:K));
end

function x = rows(mat)
	x = size(mat, 1);
end
function x = columns(mat)
	x = size(mat, 2);
end
function c = ITER_MAX()
	c = 100;
end
function c = EPSILON()
	c = 1e-4;
end

function Dmat = mkDist(Xmat, Pmat)
	Xperm = repmat(permute(Xmat, [ 3, 2, 1 ]), [ columns(Pmat), 1, 1 ]);
	Pperm = repmat(permute(Pmat, [ 2, 3, 1 ]), [ 1, columns(Xmat), 1 ]);
	Dmat = sum((Xperm - Pperm) .^2, 3);
end

function Hmat = mkAssign(Dmat)
	Hmat = zeros(size(Dmat'));
	[vals, idx] = min(Dmat);
	for i = 1:columns(Dmat)
		Hmat(i, idx(i)) = 1;
	end
end

function [L, Lvec, i] = em (Xmat, K)
	[d, N] = size(Xmat);

	% initialize Umat, Dmat, Hmat with k-means
	[Error, Evec, Iter, Retries, Umat] = kmeans(Xmat, K, @RandChoiceProto);
	Dmat = mkDist(Xmat, Umat);
	Hmat = mkAssign(Dmat);
	sigmaSquare = sum(sum(Dmat' .* Hmat)) / (N*d);

	% P(j)
	Pmat = sum(Hmat) ./ N;

	% likelihood vector
	Lvec = zeros(1, ITER_MAX);
	Lvec(1) = inf;

	i = 1;
	while ((i == 1) || abs(Lvec(i) - Lvec(i-1)) > EPSILON) && i < ITER_MAX
        i = i+1;

		% p(x|j)
		pmat = (2 * pi * sigmaSquare) ^ -(d/2) * exp( (-1 / (2*sigmaSquare)) * Dmat' );

		% assignment matrix h(i,j)
		Hmat = 	(repmat(Pmat, rows(pmat), 1) .* pmat) ./ ... % p(x|j)*P(j)	
				repmat( ...
					sum( repmat(Pmat, rows(pmat), 1) .* pmat, 2 ), 1, columns(pmat) ...
				);

		% P(j)
		Pmat = sum(Hmat) ./ N;

		% µ(j)
		Umat = Xmat * Hmat ./ repmat(sum(Hmat), rows(Xmat), 1);

		% x_i - mu_j
		Dmat = mkDist(Xmat, Umat);

		sigmaSquare = sum(sum(Hmat .* Dmat')) / (N*d);

		Lvec(i) = sum( ...
					log( ...
						sum(repmat(Pmat, rows(pmat), 1) .* pmat, 2)	... % p(x_i)
					) ...
				);
		L = Lvec(i);
	end

end



function sheet2_ex1a(Xmat)
	[L, Lvec, Iter] = em(Xmat, 4);
	plot([2:Iter], Lvec(2:Iter))
end

function sheet2_ex1b(Xmat, Kmax)
	Lvec = zeros(1, Kmax);
	for K = 2:Kmax
		[L, Lv, Iter] = em(Xmat, K);
		Lvec(K) = L;
		L
	end
	plot([2 : Kmax], Lvec(2 : end));
end


function sheet2_ex2(Xmat, Kmax)
	[d, N] = size(Xmat);
	BICvec = zeros(1, Kmax);
	for K = 2:Kmax
		[L, Lvec, Iter] = em(Xmat, K);
		M = K*d + (K-1) + 1;
		BICvec(K) = 2 * L - M * log(N);
	end
	plot([2 : Kmax], BICvec(2 : end));
end
	

