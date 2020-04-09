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
