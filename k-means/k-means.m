% turn off warnings about automatic broadcasting
warning('off', 'Octave:broadcast')


function DataMat = GenDataMat()
	dev = 0.25;
	DataMat = [ randn(2, 100).*dev,  \
		randn(2, 100).*dev + [0; 1], \
		randn(2, 100).*dev + [1; 0], \
		randn(2, 100).*dev + [1; 1] ];
end

function ProtoMat = RandProto(K, DataMat)
	ProtoMat = rand(rows(DataMat), K);
end

function ProtoMat = RandChoiceProto(K, DataMat)
	ProtoMat = DataMat(:,randperm(columns(DataMat), K));
end

function [Error, Iter, ProtoMat, Retries] = kmeans (DataMat, k, ProtoFun)
	EPSILON = 1e-4;
	ITER_MAX = 100;
	E = zeros(ITER_MAX);
	Iter = 0;
	Retries = 0;

	ProtoMat = ProtoFun(k, DataMat);

	while ++Iter < ITER_MAX && ((Iter <= 2) || abs(E(Iter-2) - E(Iter-1)) > EPSILON)

		% distance matrix
		for i = 1:columns(DataMat)
			for j = 1:columns(ProtoMat)
				DistMat(j,i) = sumsq(DataMat(:,i) - ProtoMat(:,j));
			end
		end

		AssignMat = zeros(columns(DataMat), columns(ProtoMat));
		for i = 1:columns(DistMat)
			[val, idx] = min(DistMat(:,i));
			AssignMat(i, idx) = 1;
		end

		% check if there is one prototype without assigned data vectors
		if !all(sum(AssignMat))
			ProtoMat = ProtoFun(k, DataMat);
			Iter = 0;
			Retries += 1;
			continue
		end


		for i = 1:columns(ProtoMat)
			ProtoMat(:,i) = sum(DataMat'.*AssignMat(:,i)) / sum(AssignMat(:,i));
		end

		E(Iter) = sum(sum(DistMat'.*AssignMat));

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
		DataMat = GenDataMat();
		[Error, Iterations, ProtoMat, Retries] = kmeans(DataMat, K, @RandProto);
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
		DataMat = GenDataMat();
		[Error, Iterations, ProtoMat, Retries] = kmeans(DataMat, K, @RandChoiceProto);
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
			[Error, Iterations, ProtoMat] = kmeans(DataMat, K, @RandChoiceProto);
			E(K) = min(E(K), Error);
		end
	end
	plot(E(2:end))
end

ex_a(20)
ex_b(20)
