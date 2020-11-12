# Turn a sequence of tokens (often qgrams) to a count dict for them, 
# i.e. map each word/token to the number of times it has been seen.
function countdict(qgrams)
    d = Dict{eltype(qgrams), Int32}()
    numtokens = 0
    for qg in qgrams
        numtokens += 1
        index = Base.ht_keyindex2!(d, qg)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = qg
			@inbounds d.vals[index] = d.vals[index][1] + 1
		else
			@inbounds Base._setindex!(d, 1, qg, -index)
		end
    end
    return d, numtokens
end

# We will pre-process inputs to dictionaries that keeps track of the
# words, their multiplicity (count), as well as the total count of
# word instances in the input.
# A dictionary is either of fixed or varying word length.

abstract type WordLength end
wordlength(wl::WordLength) = wordlength(typeof(wl))

struct VaryingWordLength <: WordLength end
wordlength(::Type{VaryingWordLength}) = :varying

struct FixedWordLength{Q} <: WordLength end
wordlength(::Type{FixedWordLength{Q}}) where Q = Q

abstract type AbstractWordDictionary{WL<:WordLength,K} end

# Default is that sub-types saves the key data directly in fields:
wordcounts(d::AbstractWordDictionary) = d.wordcounts
numtokens(d::AbstractWordDictionary) = d.numtokens # sum of the counts for all words/qgrams
numwords(d::AbstractWordDictionary) = d.numwords # number of unique words/qgrams

Base.length(d::AbstractWordDictionary{<:FixedWordLength,K}) where K =
    numtokens(d) > 0 ? (numtokens(d) + wordlength(d) - 1) : error("unknown length of original string")

"""
	wordlength(d::AbstractWordDictionary)

Return the length of the words in the dictionary. This is either
a positive integer (for fixed length, q-gram dictionaries) or
the symbol :varying for dictionaries where the word length varies.
"""
wordlength(d::AbstractWordDictionary{WL,K}) where {WL,K} = wordlength(WL)

wordlengthtype(qitr::QGramIterator) = FixedWordLength{qitr.q}

"""
    WordDictionary(s, q::Integer = 2)
    WordDictionary(s, itr::WordIterator = LempelZivIterator())

Creates a WordDictionary that pre-calculates (pre-counts) the words/qgrams
of a string or stream. This enables faster calculation of Dictionary/QGram 
distances.

Note that the wordlength of the dictionary must be the same as what is 
expected by the distance.

## Examples
```julia
str1, str2 = "my string", "another string"
qd1 = WordDictionary(str1, 2)
qd2 = WordDictionary(str2, 2)
evaluate(Overlap(2), qd1, qd2)
```
"""
struct WordDictionary{WL,K} <: AbstractWordDictionary{WL,K}
    wordcounts::Dict{K,Int}
    numwords::Int
    numtokens::Int
end
function WordDictionary(s::Union{AbstractString, AbstractVector}, q::Integer = 2)
    @assert q >= 1
    WordDictionary(s, qgrams(s, q))
end
function WordDictionary(s::Union{AbstractString, AbstractVector}, itr::AbstractTokenIterator)
    counts, numtokens = countdict(itr)
    WordDictionary{wordlengthtype(itr), eltype(itr)}(counts, length(counts), numtokens)
end
WordDictionary(s, qoritr::Union{Integer, AbstractTokenIterator}) = WordDictionary(collect(s), qoritr)

"""
	WordSortedVector(s, q::Integer = 2)

Creates a WordSortedVector that pre-calculates (pre-counts) the 
words/qgrams of a string or stream. This enables faster calculation of
Dictionary/QGram distances.

Since words are sorted in lexicographic order QGram distances can be 
calculated even faster than when using a WordDictionary. However, the 
sorting means that updating the counts after creation is less 
efficient. However, for most use cases WordSortedVector is preferred
over a WordDictionary.

Note that the word length must correspond with the word length of the 
distance.

## Examples
```julia
str1, str2 = "my string", "another string"
qs1 = WordSortedVector(str1, 2)
qs2 = WordSortedVector(str2, 2)
evaluate(Jaccard(2), qs1, qs2)
```
"""
struct WordSortedVector{WL,K} <: AbstractWordDictionary{WL,K}
    wordcounts::Vector{Pair{K,Int}}
    numwords::Int
    numtokens::Int
end
function WordSortedVector(s::Union{AbstractString, AbstractVector}, q::Integer = 2)
    @assert q >= 1
    WordSortedVector(s, qgrams(s, q))
end
function WordSortedVector(s::Union{AbstractString, AbstractVector}, itr::AbstractTokenIterator)
    wd = WordDictionary(s, itr)
    countpairs = collect(wordcounts(wd))
    sort!(countpairs, by = first)
    WordSortedVector{wordlengthtype(itr), eltype(itr)}(countpairs, numwords(wd), numtokens(wd))
end
WordSortedVector(s, qoritr::Union{Integer, AbstractTokenIterator}) = WordSortedVector(collect(s), qoritr)
