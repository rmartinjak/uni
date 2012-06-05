# turn off warnings about automatic broadcasting
warning("off", "Octave:broadcast")


function [iter, ProtoMat] = kmeans (DataMat, ProtoMat)
	ITER_MAX = 100;
	EPSILON = 1e-4;

	iter = 0;
	Eold = 10000;
	Enew = 1000;

	while ++iter < ITER_MAX && abs(Eold - Enew) > EPSILON

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

		for i = 1:columns(ProtoMat)
		ProtoMat(:,i) = sum(DataMat'.*AssignMat(:,i)) / sum(AssignMat(:,i));
		end

		Eold = Enew;
		Enew = sum(sum(DistMat'.*AssignMat));
	end
endfunction

# data vectors
DataMat = [ randn(2, 100).*0.5, \
	randn(2, 100).*0.5 + [0; 1], \
	randn(2, 100).*0.5 + [1; 0], \
	randn(2, 100).*0.5 + [1; 1] ];

K = 4;

ProtoMat = rand(2, K)

[Iterations, ProtoMat] = kmeans(DataMat, ProtoMat);
ProtoMat
Iterations
