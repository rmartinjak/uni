function [Error, Evec, Iter, Retries, Pmat] = kmeans (Xmat, K, ProtoFun)
	Retries = 0;
	Pmat = ProtoFun(Xmat, K);
	Evec = zeros(1, ITER_MAX);
	Evec(1) = Inf;
	i = 1;
	while ((i == 1) || abs(Evec(i) - Evec(i-1)) > EPSILON) && i < ITER_MAX
		i = i + 1;
        % calculate distance matrix
		Dmat = mkDist(Xmat, Pmat);

		% calculate vector-prototype assignment matrix
		Hmat = mkAssign(Dmat);

		% check if there is one prototype without assigned data vectors
		if ~all(sum(Hmat))
			% reinitialize prototype matrix and "restart"
			Pmat = ProtoFun(Xmat, K);
			i = 1;
			Retries = Retries + 1;
			continue
		end

		% calculate new prototypes
		Pmat = Xmat * Hmat ./ repmat(sum(Hmat), rows(Xmat), 1);

		% calculate error
		Error = sum(mean(Dmat'.*Hmat));
		Evec(i) = Error;
		end
	Iter = i;
end


function c = ITER_MAX()
	c = 100;
end
function c = EPSILON()
	c = 1e-4;
end


function x = rows(mat)
	x = size(mat, 1);
end
function x = columns(mat)
	x = size(mat, 2);
end

function Xmat = GenXmat(dev)
	Xmat = [ randn(2, 100).*dev,  ...
		randn(2, 100).*dev + repmat([0; 1], 1, 100), ...
		randn(2, 100).*dev + repmat([1; 0], 1, 100), ...
		randn(2, 100).*dev + repmat([1; 1], 1, 100) ];
end


function Pmat = RandProto(Xmat, K)
	Pmat = rand(rows(Xmat), K);
end

function Pmat = RandChoiceProto(Xmat, K)
    Perm = randperm(columns(Xmat));
	Pmat = Xmat(:,Perm(1:K));
end

function Hmat = mkAssign(Dmat)
	Hmat = zeros(size(Dmat'));
	[vals, idx] = min(Dmat);
	for i = 1:columns(Dmat)
		Hmat(i, idx(i)) = 1;
	end
end

function Dmat = mkDist(Xmat, Pmat)
	Xperm = repmat(permute(Xmat, [ 3, 2, 1 ]), [ columns(Pmat), 1, 1 ]);
	Pperm = repmat(permute(Pmat, [ 2, 3, 1 ]), [ 1, columns(Xmat), 1 ]);
	Dmat = sum((Xperm - Pperm) .^2, 3);
end
