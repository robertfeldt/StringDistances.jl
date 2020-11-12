include("helpers/iterators.jl")    # iterating over inputs, e.g. strings and vectors
include("helpers/dictionaries.jl") # dictionaries used when pre-processing
include("helpers/counters.jl")     # counter objects used for fast implementation of distances

abstract type QGramDistance <: SemiMetric end

function (dist::QGramDistance)(s1, s2)
	((s1 === missing) | (s2 === missing)) && return missing
	counter = newcounter(dist, s1, s2)
	for (n1, n2) in _count(qgrams(s1, dist.q), qgrams(s2, dist.q))
		count!(counter, n1, n2)
	end
	calculate(dist, counter)
end

# Default is to pre-process to a WordSortedVector.
preprocess(dist::QGramDistance, s) = WordSortedVector(s, dist.q)

isaligned(dist::QGramDistance, wd::AbstractWordDictionary{FWL,<:Any}) where {FWL<:FixedWordLength} =
	dist.q == wordlength(FWL)

isaligned(dist::QGramDistance, wd::AbstractWordDictionary{<:VaryingWordLength,K}) where K =
	false # for now, later we will allow varying length dictionaries if the QGramDistance allows it

function (dist::QGramDistance)(wd1::WD, wd2::WD) where {WD<:AbstractWordDictionary}
	@assert isaligned(dist, wd1)
	@assert isaligned(dist, wd2)
	counter = newcounter(dist, wd1, wd2)
	countmatches!(counter, wordcounts(wd1), wordcounts(wd2))
    calculate(dist, counter)
end

# Fallback when a distance doesn't benefit from preprocessing it need not consider
# the objects to be compared.
newcounter(d::QGramDistance, s1, s2) = newcounter(d)

"""
	QGram(q::Int)

Creates a QGram distance.

The distance corresponds to

``||v(s1, q) - v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the number of times a q-gram appears for the string s
"""
struct QGram <: QGramDistance
	q::Int
end

mutable struct SingleCounter{T, QD<:QGramDistance} <: AbstractWordMatchCounter
	shared::T
end
newcounter(d::QGram) = SingleCounter{Int, QGram}(0)

@inline count!(c::SingleCounter{Int, QGram}, n1::Integer, n2::Integer) =
	c.shared += abs(n1 - n2)

calculate(dist::QGram, c::SingleCounter{Int, QGram}) = c.shared

"""
	Cosine(q::Int)

Creates a Cosine distance.

The distance corresponds to

`` 1 - v(s1, q).v(s2, q)  / ||v(s1, q)|| * ||v(s2, q)||``

where ``v(s, q)`` denotes the vector on the space of q-grams of length q, 
that contains the  number of times a q-gram appears for the string s
"""
struct Cosine <: QGramDistance
	q::Int
end

# PP is true iff preprocessed, false if not, and nothing if we don't care
mutable struct ThreeCounters{T, QD<:QGramDistance, PP} <: AbstractWordMatchCounter
	left::T
	right::T
	shared::T
end

# No benefit to preprocessing here since the squared sum has not been preprocessed
# so needs to be calculated.
newcounter(d::Cosine) = ThreeCounters{Int, Cosine, nothing}(0, 0, 0)

@inline function count!(c::ThreeCounters{Int, Cosine, nothing}, n1::Integer, n2::Integer)
	c.left += n1^2
	c.right += n2^2
	c.shared += n1 * n2
end

calculate(d::Cosine, c::ThreeCounters{Int, Cosine, nothing}) =
	1.0 - c.shared / (sqrt(c.left) * sqrt(c.right))

"""
	Jaccard(q::Int)

Creates a Jaccard distance.

The distance corresponds to 

``1 - |Q(s1, q) ∩ Q(s2, q)| / |Q(s1, q) ∪ Q(s2, q))|``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Jaccard <: QGramDistance
	q::Int
end

calculate(d::Jaccard, c::ThreeCounters{Int, Jaccard, PP}) where PP =
	1.0 - c.shared / (c.left + c.right - c.shared)

"""
	SorensenDice(q::Int)

Creates a SorensenDice distance.

The distance corresponds to  

``1 - 2 * |Q(s1, q) ∩ Q(s2, q)|  / (|Q(s1, q)| + |Q(s2, q))|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct SorensenDice <: QGramDistance
	q::Int
end

calculate(d::SorensenDice, c::ThreeCounters{Int, SorensenDice, PP}) where PP =
	1.0 - 2.0 * c.shared / (c.left + c.right)

"""
	Overlap(q::Int)

Creates a Overlap distance.

The distance corresponds to  

``1 - |Q(s1, q) ∩ Q(s2, q)|  / min(|Q(s1, q)|, |Q(s2, q)|)``

where ``Q(s, q)``  denotes the set of q-grams of length n for the string s
"""
struct Overlap <: QGramDistance
	q::Int
end

# All of these three distances needs the numwords which we already know
# if preprocssed inputs
const IntersectionDist = Union{Jaccard, SorensenDice, Overlap}

# No need to count the words if we have preprocessed.
newcounter(d::IntersectionDist, wd1::AbstractWordDictionary, wd2::AbstractWordDictionary) = 
	ThreeCounters{Int, typeof(d), true}(numwords(wd1), numwords(wd2), 0)

newcounter(d::IntersectionDist, s1, s2) = ThreeCounters{Int, typeof(d), false}(0, 0, 0)

@inline function count!(c::ThreeCounters{Int, QD, PP}, n1::Integer, n2::Integer) where {PP,QD<:IntersectionDist}
	if !PP
		c.left += (n1 > 0)
		c.right += (n2 > 0)
	end
	c.shared += (n1 > 0) & (n2 > 0)
end

calculate(d::Overlap, c::ThreeCounters{Int, Overlap, PP}) where PP =
	1.0 - c.shared / min(c.left, c.right)

"""
	MorisitaOverlap(q::Int)

Creates a MorisitaOverlap distance, a general, statistical measure of
dispersion which can also be used on dictionaries such as created
from q-grams. See https://en.wikipedia.org/wiki/Morisita%27s_overlap_index
This is more fine-grained than many of the other QGramDistances since
it is based on the counts per q-gram rather than only which q-grams are
in the strings.

The distance corresponds to

``(2 * sum(m(s1) .* m(s2)) / (sum(m(s1).^2)*M(s2)/M(s1) + sum(m(s2).^2)*M(s1)/M(s2))``

where ``m(s)`` is the vector of q-gram counts for string ``s`` and ``M(s)`` is the
sum of those counts.
"""
struct MorisitaOverlap <: QGramDistance
	q::Int
end

mutable struct FiveCounters{T, QD<:QGramDistance, PP} <: AbstractWordMatchCounter
	leftsum::T    # sum(m(s1)), i.e. numtokens for s1 (if preprocessed)
	rightsum::T   # sum(m(s2)), i.e. numtokens for s2 (if preprocessed)
	leftsq::T     # sum(m(s1).^2)
	rightsq::T    # sum(m(s2).^2)
	shared::T     # sum(m(s1) .* m(s2))
end

newcounter(d::MorisitaOverlap, s1, s2) = FiveCounters{Int, MorisitaOverlap, false}(0, 0, 0, 0, 0)

newcounter(d::MorisitaOverlap, wd1::AbstractWordDictionary, wd2::AbstractWordDictionary) = 
	FiveCounters{Int, MorisitaOverlap, true}(numtokens(wd1), numtokens(wd2), 0, 0, 0)

# To consider: For WordSortedVector we actually have the multiplicity vectors themselves
# so we could calculate also leftsq and rightsq in newcounter. But complicates the code 
# so skipped for now.

@inline function count!(c::FiveCounters{Int, MorisitaOverlap, PP}, n1::Integer, n2::Integer) where PP
	if !PP
		c.leftsum += n1
		c.rightsum += n2
	end
	c.leftsq += (n1^2)
	c.rightsq += (n2^2)
	c.shared += (n1 * n2)
end

calculate(d::MorisitaOverlap, c::FiveCounters{Int, MorisitaOverlap, PP}) where PP =
	1.0 - ((2 * c.shared) / (c.leftsq*c.rightsum/c.leftsum + c.rightsq*c.leftsum/c.rightsum))

"""
	NMD(q::Int)

Creates a NMD (Normalized Multiset Distance) as introduced by Besiris and
Zigouris 2013. The goal with this distance is to behave similarly to a normalized
compression distance without having to do any actual compression (and thus being
faster to compute).

The distance corresponds to

``(sum(max.(m(s1), m(s2)) - min(M(s1), M(s2))) / max(M(s1), M(s2))``

where ``m(s)`` is the vector of q-gram counts for string ``s`` and ``M(s)`` is the
sum of those counts.

For details see:
https://www.sciencedirect.com/science/article/pii/S1047320313001417
"""
struct NMD <: QGramDistance
	q::Int
end

# Default is to count left, right, and shared but if we have preprocessed 
# we only need to count shared since the number of tokens is known.
newcounter(d::NMD, s1, s2) = ThreeCounters{Int, NMD, false}(0, 0, 0)
newcounter(d::NMD, wd1::AbstractWordDictionary, wd2::AbstractWordDictionary) = 
	ThreeCounters{Int, NMD, true}(numtokens(wd1), numtokens(wd2), 0)

@inline function count!(c::ThreeCounters{Int, NMD, PP}, n1::Integer, n2::Integer) where PP
	if !PP
		c.left += n1
		c.right += n2
	end
	c.shared += max(n1, n2)
end

calculate(d::NMD, c::ThreeCounters{Int, NMD, PP}) where PP =
	(c.shared - min(c.left, c.right)) / max(c.left, c.right)
