source ab2_kmeans.m

function em (Xmat, K)
	N = columns(Xmat);
	d = rows(Xmat);

	[Error, Iter, Umat, Retries] = kmeans(Xmat, 4, @RandChoiceProto);

	Dmat = mkDist(Xmat, Umat);
	Hmat = mkAssign(Dmat);

	LastError = Inf;
	Error = Inf;
	Iter = 0;
	while ++Iter < ITER_MAX && ((Iter <= 2) || abs(LastError - Error) > EPSILON)


		Pmat = sum(Hmat) ./ N;
		Umat = Xmat * Hmat ./ sum(Hmat);
		sigmasquare = sum(sum(

	end
end
