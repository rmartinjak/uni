# turn off warnings about automatic broadcasting
warning("off", "Octave:broadcast")


function [E, iter, ProtoMat] = kmeans (DataMat, ProtoMat)
	ITER_MAX = 100;
	EPSILON = 1e-4;

	iter = 0;

	while ++iter < ITER_MAX && ((iter == 1) || abs(E(iter) - E(iter-1)) > EPSILON)

		# distance matrix
		for i = 1:columns(DataMat)
		for k = 1:columns(ProtoMat)
		DistMat(k,i) = sumsq(DataMat(:,i) - ProtoMat(:,k));
		end
		end

		AssignMat = zeros(columns(DataMat), columns(ProtoMat));

		for i = 1:columns(DistMat)
		[val, idx] = min(DistMat(:,i));
		AssignMat(i, idx) = 1;
		end

		if !all(sum(AssignMat))
			1+1
		endif


		for i = 1:columns(ProtoMat)
		ProtoMat(:,i) = sum(DataMat'.*AssignMat(:,i)) / sum(AssignMat(:,i));
		end

		E(iter) = sum(sum(DistMat'.*AssignMat));
	end
endfunction

# data vectors
DataMat = [ randn(2, 100).*0.5, \
	randn(2, 100).*0.5 + [0; 1], \
	randn(2, 100).*0.5 + [1; 0], \
	randn(2, 100).*0.5 + [1; 1] ];

K = 4;

ProtoMat = rand(2, K)

[Error, Iterations, ProtoMat] = kmeans(DataMat, ProtoMat);
ProtoMat
Iterations
Error


# for K = 2:10
# ProtoMat = randperm(columns(DataMat), K);
# [Error, Iterations, ProtoMat] = kmeans(DataMat, ProtoMat);
# ProtoMat
# Iterations
# Error
