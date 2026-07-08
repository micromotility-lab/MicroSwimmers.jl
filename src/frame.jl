    struct Frame{T}
        location::SVector{3,T}
        orientation::SMatrix{3,3,T}
    end

    Frame{T}() where {T} = Frame(zero(SVector{3,T}), SMatrix{3,3,T}(I))
    Frame(l, o) = Frame(SVector{3,eltype(l)}(l), SMatrix{3,3,eltype(o)}(o))

    # SE(3) composition: (parent ∘ child) gives child's pose in world
    @inline Base.:*(P::Frame, C::Frame) =
        Frame(P.location + P.orientation * C.location, P.orientation * C.orientation)

    # act on a body-frame point
    @inline (F::Frame)(X::SVector{3}) = F.location + F.orientation * X

    @inline Base.inv(F::Frame) =
    Frame(-F.orientation' * F.location, F.orientation')




