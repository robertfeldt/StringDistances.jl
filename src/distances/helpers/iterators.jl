abstract type AbstractTokenIterator{S <: Union{AbstractString, AbstractVector}} end

struct QGramIterator{S <: Union{AbstractString, AbstractVector}} <: AbstractTokenIterator{S}
	s::S   # Collection
	q::Int # Length of Qgram
end
Base.length(qgram::QGramIterator) = max(length(qgram.s) - qgram.q + 1, 0)

# q-grams of AbstractString
function Base.iterate(qgram::QGramIterator{<: AbstractString}, 
	state = (1, nextind(qgram.s, 0, qgram.q)))
	istart, iend = state
	iend > ncodeunits(qgram.s) && return nothing
	element = SubString(qgram.s, istart, iend)
	nextstate = nextind(qgram.s, istart), nextind(qgram.s, iend)
	element, nextstate
end
Base.eltype(qgram::QGramIterator{SubString{S}}) where {S} = SubString{S}
Base.eltype(qgram::QGramIterator{S}) where {S <: AbstractString} = SubString{S}
qgrams(s::AbstractString, q::Integer) = QGramIterator(s, q)

#q-grams of AbstractVector
# Alternatively, I could also use partition in IterTools but it creates a vector for each iteration
# so it does not seem to be worth it.
function Base.iterate(qgram::QGramIterator{<: AbstractVector}, state = firstindex(qgram.s))
	state + qgram.q - 1 > lastindex(qgram.s) && return nothing
	view(qgram.s, state:(state + qgram.q - 1)), state + 1
end
Base.eltype(qgram::QGramIterator{<: AbstractVector}) = typeof(first(qgram))
qgrams(s::AbstractVector, q::Integer) = QGramIterator(s, q)
qgrams(s, q::Integer) = QGramIterator(collect(s), q)

@doc """
Return an iterator corresponding to the the q-gram of an iterator. 
When the iterator is a String, qgrams are SubStrings.

### Arguments
* `s` iterator
* `q::Integer`: length of q-gram

## Examples
```julia
for x in qgrams("hello", 2)
	println(x)
end
```
""" 
qgrams

# For two iterators s1 and s2, that define a length and eltype method,
# this returns an iterator that,
# for each element in s1 ∪ s2, returns (numbers of times it appears in s1, numbers of times it appears in s2)
function _count(s1, s2)
	K = promote_type(eltype(s1), eltype(s2))
	d = Dict{K, Tuple{Int, Int}}()
	sizehint!(d, length(s1) + length(s2))
	# I use a faster way to change a dictionary key
	# see setindex! in https://github.com/JuliaLang/julia/blob/master/base/dict.jl#L380
	for x1 in s1
		index = Base.ht_keyindex2!(d, x1)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = x1
			@inbounds d.vals[index] = (d.vals[index][1] + 1, 0)
		else
			@inbounds Base._setindex!(d, (1, 0), x1, -index)
		end
	end
	for x2 in s2
		index = Base.ht_keyindex2!(d, x2)
		if index > 0
			d.age += 1
			@inbounds d.keys[index] = x2
			@inbounds d.vals[index] = (d.vals[index][1], d.vals[index][2] + 1)
		else
			@inbounds Base._setindex!(d, (0, 1), x2, -index)
		end
	end
	return values(d)
end
