using StringDistances, Unicode, Test, Random

@testset "Word dictionaries" begin

@testset "length" begin
    @test length(WordDictionary("arn", 1)) == 3
    @test length(WordDictionary("arn", 2)) == 3
    @test length(WordDictionary("arn", 3)) == 3
    @test_throws ErrorException length(WordDictionary("arn", 4))

    for _ in 1:100
        s = randstring(rand(10:100))
        for WD in [WordDictionary, WordSortedVector]
            wd1 = WD(s, 1)
            @test length(wd1) == length(s)

            wd2 = WD(s, 2)
            @test length(wd2) == length(s)
        end
    end
end

@testset "WordDictionary and WordSortedVector counts qgrams" begin
    # To get something we can more easily compare to:
    stringify(p::Pair{<:AbstractString, <:Integer}) = (string(first(p)), last(p))
    stringify(p::Pair{V, <:Integer}) where {S<:AbstractString,V<:AbstractVector{S}} = (map(string, first(p)), last(p))
    sortedcounts(qc) = sort(collect(StringDistances.wordcounts(qc)), by = first)
    totuples(qc) = map(stringify, sortedcounts(qc))

    s1, s2   = "arnearne", "arnebeda"

    qd1, qd2 = WordDictionary(s1, 2), WordDictionary(s2, 2)
    @test totuples(qd1) == [("ar", 2), ("ea", 1), ("ne", 2), ("rn", 2)]
    @test totuples(qd2) == [("ar", 1), ("be", 1), ("da", 1), ("eb", 1), ("ed", 1), ("ne", 1), ("rn", 1)]

    qc1, qc2 = WordSortedVector(s1, 2), WordSortedVector(s2, 2)
    @test totuples(qc1) == [("ar", 2), ("ea", 1), ("ne", 2), ("rn", 2)]
    @test totuples(qc2) == [("ar", 1), ("be", 1), ("da", 1), ("eb", 1), ("ed", 1), ("ne", 1), ("rn", 1)]

    s3 = "rgówów"
    qd3a = WordDictionary(s3, 2)
    @test totuples(qd3a) == [("gó", 1), ("rg", 1), ("wó", 1), ("ów", 2)]

    qd3b = WordDictionary(graphemes(s3), 2)
    @test totuples(qd3b) == [(["g", "ó"], 1), (["r", "g"], 1), (["w", "ó"], 1), (["ó", "w"], 2)]

    qc3a = WordSortedVector(s3, 2)
    @test totuples(qc3a) == [("gó", 1), ("rg", 1), ("wó", 1), ("ów", 2)]

    qd3b = WordDictionary(graphemes(s3), 2)
    @test totuples(qd3b) == [(["g", "ó"], 1), (["r", "g"], 1), (["w", "ó"], 1), (["ó", "w"], 2)]
end

end