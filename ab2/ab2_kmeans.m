
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

function Xmat = GenXmat(dev)
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
	Dmat = sumsq((Xperm - Pperm), 3);
end

function [Error, Evec, Iter, Retries, Pmat] = kmeans (Xmat, K, ProtoFun)
	Retries = 0;
	Pmat = ProtoFun(Xmat, K);
	Evec = zeros(1, ITER_MAX);
	Evec(1) = Inf;
	i = 1;
	while ((i == 1) || abs(Evec(i) - Evec(i-1)) > EPSILON) && ++i < ITER_MAX
		% calculate distance matrix
		Dmat = mkDist(Xmat, Pmat);

		% calculate vector-prototype assignment matrix
		Hmat = mkAssign(Dmat);

		% check if there is one prototype without assigned data vectors
		if !all(sum(Hmat))
			% reinitialize prototype matrix and "restart"
			Pmat = ProtoFun(Xmat, K);
			i = 1;
			Retries++;
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


% 1a) / 1b)

function [Error, Retries] = sheet1_ex_a_b (K, Iterations, ProtoFun)
	Evec = zeros(1, Iterations);
	Retries = 0;
	for i = 1:Iterations
		Xmat = GenXmat(0.25);
		[Error, _, _, Ret, _] = kmeans (Xmat, K, ProtoFun);
		Evec(i) = Error;
		Retries += Ret;
	end
	Error = mean(Evec);
end

function [Error, Retries] = sheet1_ex_a (K, Iterations)
	% plot quadratic error
	Xmat = GenXmat(0.25);
	[_, Evec, Iter, _, _] = kmeans(Xmat, K, @RandProto);
	plot([2:Iter], Evec(2:Iter))

	% expectation: dev^2 = 0.0625
	[Error, Retries] = sheet1_ex_a_b(K, Iterations, @RandProto)
end
function [Error, Retries] = sheet1_ex_b (K, Iterations)
	[Error, Retries] = sheet1_ex_a_b(K, Iterations, @RandChoiceProto)
end


% 1c) / 2

function sheet1_ex_c (Xmat, Kmax, Iterations)
	Evec = zeros(1, Kmax);
	for K = 2:Kmax
		Evec(K) = Inf;
		for i = 1:Iterations
			[Error, _, _] = kmeans(Xmat, K, @RandChoiceProto);
			Evec(K) = min(Evec(K), Error);
		end
	end
	plot([2 : Kmax], Evec(2 : end));
end


%############%
% Exercise 1 %
%############%
% sheet1_ex_a(4, 50);
% sheet1_ex_b(4, 50);
% sheet1_ex_c(GenXmat(0.25), 10, 20);

%############%
% Exercise 2 %
%############%
% load GexprData1.mat;
% sheet1_ex_c(Xmat, 10, 20);
