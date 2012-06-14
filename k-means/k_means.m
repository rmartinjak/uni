% turn off warnings about automatic broadcasting
warning('off', 'Octave:broadcast')


function Xmat = GenXmat()
	dev = 0.25;
	Xmat = [ randn(2, 100).*dev,  \
		randn(2, 100).*dev + repmat([0; 1], 1, 100), \
		randn(2, 100).*dev + repmat([1; 0], 1, 100), \
		randn(2, 100).*dev + repmat([1; 1], 1, 100) ];
end

function Pmat = RandProto(Xmat, K)
	Pmat = rand(rows(Xmat), K);
end

function Pmat = RandChoiceProto(Xmat, K)
	Pmat = Xmat(:,randperm(columns(Xmat))(1:K));
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
	Dmat = sum((Xperm - Pperm) .^ 2, 3);
end

function [Error, Iter, Pmat, Retries] = kmeans (Xmat, K, ProtoFun)
	EPSILON = 1e-4;
	ITER_MAX = 100;
	E = zeros(ITER_MAX);
	Iter = 0;
	Retries = 0;

	Pmat = ProtoFun(Xmat, K);

	while ++Iter < ITER_MAX && ((Iter <= 2) || abs(E(Iter-2) - E(Iter-1)) > EPSILON)

		Dmat = mkDist(Xmat, Pmat);
		Hmat = mkAssign(Dmat);

		% check if there is one prototype without assigned data vectors
		if !all(sum(Hmat))
			Pmat = ProtoFun(Xmat, K);
			Iter = 0;
			Retries += 1;
			continue
		end

		% calculate new prototypes
		Pmat = Xmat * Hmat ./ repmat(sum(Hmat), rows(Xmat), 1);

		E(Iter) = sum(mean(Dmat'.*Hmat));

		end

	Error = E(Iter-1);
end


%===========
% 1a) / 1b)
%===========

function [Error, Retries] = ex_a_b (K, iterations, RandProto)
	E = zeros(1, iterations);
	Retries = 0;
	for Iter = 1:iterations
		Xmat = GenXmat();
		[Err, Iterations, Pmat, Ret] = kmeans(Xmat, K, RandProto);
		E(Iter) = Err;
		Retries += Ret;
	end
	Error = mean(E);
end

function [Error, Retries] = ex_a (K, iterations)
	[Error, Retries] = ex_a_b(K, iterations, @RandProto);
end
function [Error, Retries] = ex_b (K, iterations)
	[Error, Retries] = ex_a_b(K, iterations, @RandChoiceProto);
end


%=========
% 1c) / 2
%=========

function ex_c (Xmat, iterations)
	E = zeros(1, 10);
	for K = 2:columns(E)
		E(K) = Inf;
		for i = 1:iterations
			[Error, Iterations, Pmat] = kmeans(Xmat, K, @RandChoiceProto);
			E(K) = min(E(K), Error);
		end
	end
	plot([2: columns(E)], E(2:end));
end
