% turn off warnings about automatic broadcasting
warning("off", "Octave:broadcast")

function ProtoMat = RandProto(K)
	ProtoMat = rand(2, K);
endfunction

function ProtoMat = RandChoiceProto(K, DataMat)
	ProtoMat = DataMat(:,randperm(columns(DataMat), K));
endfunction

function [Error, iter, ProtoMat] = kmeans (DataMat, k, ProtoFun)
	ITER_MAX = 100;
	EPSILON = 1e-4;

	E = zeros(ITER_MAX);
	iter = 0;

	ProtoMat = ProtoFun(k, DataMat);

	while ++iter < ITER_MAX && ((iter <= 2) || abs(E(iter-2) - E(iter-1)) > EPSILON)

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
			printf("zomg!\n");
			ProtoMat = ProtoFun(k, DataMat);
			iter = 0;
			continue
		endif


		for i = 1:columns(ProtoMat)
			ProtoMat(:,i) = \
			sum(DataMat'.*AssignMat(:,i)) / sum(AssignMat(:,i));
		end

		 E(iter) = sum(sum(DistMat'.*AssignMat));

	end

	Error = E(iter-1);
endfunction

% data vectors
DataMat = [ randn(2, 100).*0.5, \
	randn(2, 100).*0.5 + [0; 1], \
	randn(2, 100).*0.5 + [1; 0], \
	randn(2, 100).*0.5 + [1; 1] ];

K = 4;

[Error, Iterations, ProtoMat] = kmeans(DataMat, K, @RandProto);
ProtoMat
Iterations
%plot(Error)


E = zeros(1, 10);
for K = 2:columns(E)
	E(K) = Inf;
	for i = 1:20
		[Error, Iterations, ProtoMat] = kmeans(DataMat, K, @RandChoiceProto);
		E(K) = min(E(K), Error);
	end
	E
end
plot(E(2:end))
