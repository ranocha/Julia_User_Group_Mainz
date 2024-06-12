using BenchmarkTools
function diff1!(Y,X)
    for i in eachindex(Y)
        Y[i] = X[i+1]-X[i]
    end
end
function diff2!(Y,X)
    Y .= X[2:end] .- X[1:end-1];
end
function diff3!(Y,X)
    @views Y .= X[2:end] .- X[1:end-1];
end
function diff4!(Y,X)
    @views for i in eachindex(Y)
        Y[i] = X[i+1] - X[i]
    end
end
diff5(X) = X[2:end] .- X[1:end-1];
@views diff6(X) = X[2:end] .- X[1:end-1];
#Define variables and Use functions
X = [1:10.0;];
Y = zeros(length(X)-1);
#Run btime to check allocations of various diff functions
@btime diff1!(Y,X);
@btime diff2!(Y,X);
@btime diff3!(Y,X);
@btime diff4!(Y,X);
@btime Y = diff5(X);
@btime Y.= diff6(X);
@btime Y = diff6(X);
@btime Y = diff(X);
@btime Y.= diff(X);
