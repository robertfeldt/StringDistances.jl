using Pkg, BenchmarkTools, Random

if length(ARGS) > 2 && ARGS[3] == "local"
    Pkg.activate(pwd())
end
using StringDistances

println("isdefined(StringDistances, :QGramDict) = ", isdefined(StringDistances, :QGramDict))
println("isdefined(StringDistances, :WordDictionary) = ", isdefined(StringDistances, :WordDictionary))

N = if length(ARGS) > 0
    try
        parse(Int, ARGS[1])
    catch _
        100
    end
else
    100 # default value
end

Maxlength = if length(ARGS) > 1
    try
        parse(Int, ARGS[2])
    catch _
        100
    end
else
    100 # default value
end

# If there are strings already cached to disk we start with them and only
# add new ones if needed.
using Serialization
const CacheFile = joinpath(@__DIR__(), "perfteststrings_$(Maxlength).juliabin")
SaveCache = false

S = if isfile(CacheFile)
    try
        res = deserialize(CacheFile)
        println("Read $(length(res)) strings from cache file: $CacheFile")
        res
    catch err
        String[]
    end
else
    println("Creating $N random strings.")
    SaveCache = true
    String[randstring(rand(3:Maxlength)) for _ in 1:N]
    SaveCache = true
    String[randstring(rand(3:Maxlength)) for _ in 1:N]
end

if length(S) < N
    for i in (length(S)+1):N
        push!(S, randstring(rand(3:Maxlength)))
    end
    SaveCache = true
end

if SaveCache
    println("Saving cache file with $(length(S)) strings: $CacheFile")
    serialize(CacheFile, S)
end

S = S[1:N]

DT = if length(ARGS) == 4
    StringDistances.eval(Symbol(ARGS[4]))
else
    Cosine
end

dist = DT(2)
println("Distance = ", dist)
println("For ", Threads.nthreads(), " threads and ", length(S), " strings of max length ", Maxlength, ":")
t1 = @belapsed dm1 = pairwise(dist, S; preprocess = false)
t2 = @belapsed dm2 = pairwise(dist, S; preprocess = true)

println("  - time WITHOUT pre-calculation: ", round(t1, digits = 3))
println("  - time WITH    pre-calculation: ", round(t2, digits = 3))
println("  - speedup with pre-calculation: ", round(t1/t2, digits = 1))
