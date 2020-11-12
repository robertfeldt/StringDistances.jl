# To implement the distances we will count word/qgram matches
# between strings or pre-calculated AbstractWordMatchCounter objects.
abstract type AbstractWordMatchCounter end
@inline count!(c::AbstractWordMatchCounter, qg, n1::Integer, n2::Integer) = count!(c, n1, n2)
@inline count!(c::AbstractWordMatchCounter, n1::Integer, n2::Integer) = nothing

# The fastest version uses pairs of words and their integer counts.
function countmatches!(mc::AbstractWordMatchCounter, d1::Vector{Pair{K,I}}, d2::Vector{Pair{K,I}}) where {K,I<:Integer}
    i1 = i2 = 1
    while i1 <= length(d1) || i2 <= length(d2)
        if i2 > length(d2)
			for i in i1:length(d1)
				@inbounds count!(mc, d1[i][1], d1[i][2], 0)
            end
            return
        elseif i1 > length(d1)
			for i in i2:length(d2)
				@inbounds count!(mc, d2[i][1], 0, d2[i][2])
            end
            return
        end
        @inbounds k1, n1 = d1[i1]
        @inbounds k2, n2 = d2[i2]
        cmpval = Base.cmp(k1, k2)
		if cmpval == -1 # k1 < k2
			count!(mc, k1, n1, 0)
            i1 += 1
        elseif cmpval == +1 # k2 < k1
			count!(mc, k2, 0, n2)
            i2 += 1
		else
			count!(mc, k1, n1, n2)
            i1 += 1
            i2 += 1
        end
    end
end

function countmatches!(mc::AbstractWordMatchCounter, d1::Dict{K,I}, d2::Dict{K,I}) where {K,I<:Integer}
    for (k1, c1) in d1
        index = Base.ht_keyindex2!(d2, k1)
		if index > 0
			count!(mc, k1, c1, d2.vals[index])
		else
			count!(mc, k1, c1, 0)
        end
    end
    for (k2, c2) in d2
        index = Base.ht_keyindex2!(d1, k2)
		if index <= 0
			count!(mc, k2, 0, c2)
        end
    end
end
