% turn off warnings about automatic broadcasting
warning('off', 'Octave:broadcast')


function Xmat = GenXmat()
	dev = 0.25;
	Xmat = [ randn(2, 100).*dev,  \
		randn(2, 100).*dev + [0; 1], \
		randn(2, 100).*dev + [1; 0], \
		randn(2, 100).*dev + [1; 1] ];
end

function Pmat = RandProto(K, Xmat)
	Pmat = rand(rows(Xmat), K);
end

function Pmat = RandChoiceProto(K, Xmat)
	Pmat = Xmat(:,randperm(columns(Xmat), K));
end

function Hmat = mkAssign(Dmat)
	Hmat = zeros(size(Dmat'));
	for i = 1:columns(Dmat)
		[val, idx] = min(Dmat(:,i));
		Hmat(i, idx) = 1;
	end
end

function Dmat = mkDist(Xmat, Pmat)
	for i = 1:columns(Xmat)
		for j = 1:columns(Pmat)
			Dmat(j,i) = sumsq(Xmat(:,i) - Pmat(:,j));
		end
	end
end

function [Error, Iter, Pmat, Retries] = kmeans (Xmat, k, ProtoFun)
	EPSILON = 1e-4;
	ITER_MAX = 100;
	E = zeros(ITER_MAX);
	Iter = 0;
	Retries = 0;

	Pmat = ProtoFun(k, Xmat);

	while ++Iter < ITER_MAX && ((Iter <= 2) || abs(E(Iter-2) - E(Iter-1)) > EPSILON)

		Dmat = mkDist(Xmat, Pmat);
		Hmat = mkAssign(Dmat);

		% check if there is one prototype without assigned data vectors
		if !all(sum(Hmat))
			Pmat = ProtoFun(k, Xmat);
			Iter = 0;
			Retries += 1;
			continue
		end

		% calculate new prototypes
		Pmat = Xmat * Hmat ./ sum(Hmat);

		E(Iter) = sum(sum(Dmat'.*Hmat));

		end

	Error = E(Iter-1);
end


%====
% a)
%====

function ex_a (iterations)
	K = 4;
	E = zeros(1, iterations);
	Ret = 0;
	for Iter = 1:iterations
		Xmat = GenXmat();
		[Error, Iterations, Pmat, Retries] = kmeans(Xmat, K, @RandProto);
		E(Iter) = Error;
		Ret += Retries;
	end

	mean(E)
	Ret
end



%====
% b)
%====

function ex_b (iterations)
	K = 4;
	E = zeros(1, iterations);
	Ret = 0;
	for Iter = 1:iterations
		Xmat = GenXmat();
		[Error, Iterations, Pmat, Retries] = kmeans(Xmat, K, @RandChoiceProto);
		E(Iter) = Error;
		Ret += Retries;
	end

	mean(E)
	Ret
end


%====
% c)
%====

function ex_c (iterations)
	E = zeros(1, 10);
	for K = 2:columns(E)
		E(K) = Inf;
		for i = 1:iterations
			[Error, Iterations, Pmat] = kmeans(Xmat, K, @RandChoiceProto);
			E(K) = min(E(K), Error);
		end
	end
	plot(E(2:end))
end

