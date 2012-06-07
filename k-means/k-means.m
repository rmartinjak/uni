% turn off warnings about automatic broadcasting
warning('off', 'Octave:broadcast')

function DataMat = GenDataMat()
	DataMat = [ randn(2, 100).*0.5, \
		randn(2, 100).*0.5 + [0; 1], \
		randn(2, 100).*0.5 + [1; 0], \
		randn(2, 100).*0.5 + [1; 1] ];
end

function ProtoMat = RandProto(K)
	ProtoMat = rand(2, K);
endfunction

function ProtoMat = RandChoiceProto(K, DataMat)
	ProtoMat = DataMat(:,randperm(columns(DataMat), K));
endfunction

function [Error, Iter, ProtoMat, Retries] = kmeans (DataMat, k, ProtoFun)
	ITER_MAX = 100;
	EPSILON = 1e-4;

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
			++Retries;
			continue
		endif


		for i = 1:columns(ProtoMat)
			ProtoMat(:,i) = sum(DataMat'.*AssignMat(:,i)) / sum(AssignMat(:,i));
		end

		 E(Iter) = sum(sum(DistMat'.*AssignMat));

	end

	Error = E(Iter-1);
endfunction


%====
% a)
%====

K = 4;
E = zeros(1, 50);
Ret = 0;
for Iter = 1:50
	DataMat = GenDataMat();
	[Error, Iterations, ProtoMat, Retries] = kmeans(DataMat, K, @RandProto);
	E(Iter) = Error;
	Ret += Retries;
end

mean(E)
Ret



%====
% b)
%====

K = 4;
E = zeros(1, 50);
Ret = 0;
for Iter = 1:50
	DataMat = GenDataMat();
	[Error, Iterations, ProtoMat, Retries] = kmeans(DataMat, K, @RandChoiceProto);
	E(Iter) = Error;
	Ret += Retries;
end

mean(E)
Ret


%===
% c)
%===

E = zeros(1, 10);
for K = 2:columns(E)
	E(K) = Inf;
	for i = 1:20
		[Error, Iterations, ProtoMat] = kmeans(DataMat, K, @RandChoiceProto);
		E(K) = min(E(K), Error);
	end
end
plot(E(2:end))
