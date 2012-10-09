function ab2_kmeans
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



% 1a) / 1b)

function [Error, Retries] = sheet1_ex_a_b (K, Iterations, ProtoFun)
	Evec = zeros(1, Iterations);
	Retries = 0;
	for i = 1:Iterations
		Xmat = GenXmat(0.25);
        [Error, Evec, Iter, Ret, Pmat] = kmeans (Xmat, K, ProtoFun);
		Evec(i) = Error;
		Retries = Retries + Ret;
	end
	Error = mean(Evec);
end

function [Error, Retries] = sheet1_ex_a (K, Iterations)
	% plot quadratic error
	Xmat = GenXmat(0.25);
	[Error, Evec, Iter, Ret, Pmat] = kmeans(Xmat, K, @RandProto);
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
			[Error, Evec, Iter, Ret, Pmat]= kmeans(Xmat, K, @RandChoiceProto);
			Evec(K) = min(Evec(K), Error);
		end
	end
	plot([2 : Kmax], Evec(2 : end));
end



