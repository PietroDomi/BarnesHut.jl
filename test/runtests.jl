if ! isdir("test-gifs")
    mkdir("test-gifs")
end

include("brute-force-test.jl")
include("barnes-hut-test.jl")
println("Animations tested and available at ./test/test-gifs")

println("All tests passed!")