# just some Julia for plotting the benchmarks and examples

using Plots
using DelimitedFiles
using Statistics
A = readdlm("areaData.txt",',');
T = readdlm("runtimeData.txt",',');

D = collect(0:size(A,1))
inds = [0,8*4 .^D...].+1
Areas = [Float64.(A[i,inds[i]:inds[i+1]-1]) for i in 1:size(D,1)-1];

histogram(Areas[11],normed=true,label="",title="Histogram of Trixel areas at depth 10")
savefig("Depth10Histogram.pdf")
p1 = scatter(Int.(collect(0:size(Areas,1)-1)),mean.(Areas),ribbon=std.(Areas),label="",title="Mean and std Trixel Area at Depth d",xlabel="Depth d", ylabel="Area")
p2 = scatter(Int.(collect(0:size(Areas,1)-1)),mean.(Areas),ribbon=std.(Areas),label="",title="Mean and std Trixel Area at Depth d",xlabel="Depth d",yaxis=:log10,ylabel="Log10 Area")
plot(p1,p2,size=(900,600))
savefig("MeanStdTrixelArea.pdf")
plot(collect(0:(size(T',1)-2)).+1,T[1,1:end-1],label="",xaxis=:log10,yaxis=:log10,title = "Runtime loglog",xlabel="depth+1",ylabel="time, s")
savefig("TimeComplexity.pdf")

mutable struct Trixel{T<:AbstractArray}
    vertices::T
    id::String
end
function Trixel(;vertices=zeros(2,3))
    if size(vertices,1) != 3
        println("Must have exactly 3 vertices")
        return Triangle(zeros(2,3),"")
    end
    return Trixel(vertices,"")
end
function DrawTrixel!(T::Trixel;subplot=nothing)
    if subplot == nothing
        plot!([T.vertices[1,1],T.vertices[1,2]],
              [T.vertices[2,1],T.vertices[2,2]],
              [T.vertices[3,1],T.vertices[3,2]],label="")
        plot!([T.vertices[1,1],T.vertices[1,3]],
          [T.vertices[2,1],T.vertices[2,3]],
          [T.vertices[3,1],T.vertices[3,3]],label="")
        plot!([T.vertices[1,2],T.vertices[1,3]],
          [T.vertices[2,2],T.vertices[2,3]],
          [T.vertices[3,2],T.vertices[3,3]],label="")
    else
        plot!([T.vertices[1,1],T.vertices[1,2]],
              [T.vertices[2,1],T.vertices[2,2]],
              [T.vertices[3,1],T.vertices[3,2]],label="",subplot=subplot)
        plot!([T.vertices[1,1],T.vertices[1,3]],
          [T.vertices[2,1],T.vertices[2,3]],
          [T.vertices[3,1],T.vertices[3,3]],label="",subplot=subplot)
        plot!([T.vertices[1,2],T.vertices[1,3]],
          [T.vertices[2,2],T.vertices[2,3]],
          [T.vertices[3,2],T.vertices[3,3]],label="",subplot=subplot)
    end
end

p = plot(layout=5,aspect_ratio=:equal)
for i in 0:4
    HTM = readdlm("HTM-$i.txt",',')
    T = []
    for j in 1:size(HTM,1)
        v = Float64.(HTM[j,1:9])
        trix = Trixel(vertices=reshape(v,(3,3)))
        push!(T,trix)
        DrawTrixel!(trix,subplot=i+1)
    end
    plot!(title="HTM: depth = $i",subplot=i+1)
end
plot!(xaxis=nothing,yaxis=nothing,zaxis=nothing)
savefig(p,"ExampleHTMs.pdf")
