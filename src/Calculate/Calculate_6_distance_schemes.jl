using StaticArrays
using LinearAlgebra

export calculate_distance

function calculate_distance(point1::SVector{N, T}, point2::SVector{N, T}, method::String = "euclidean") where {N, T <: Real}
    if method == "euclidean"
        return norm(point1 - point2)
    elseif method == "manhattan"
        return sum(abs.(point1 - point2)) 
    elseif method == "chebysev"
        return maximum(abs.(point1 - point2))
    elseif method == "euclidean_squared"
        # return dot(point1 - point2, point1 - point2)

        return sum((point1 - point2).^2)
    else
        error("Unsupported method: $method.")
    end
end