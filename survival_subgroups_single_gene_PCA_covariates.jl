using Distributed
using ArgParse

@everywhere using JLD
@everywhere using StatsBase
@everywhere using SharedArrays



function userArgs()
    arglist = ArgParseSettings()
    @add_arg_table arglist begin
        "--cpus","-c"
            help = "CPU number for parallel computation"
            arg_type = Int
            default = 8
        "--iteration","-n"
            help = "Iteration times"
            arg_type = Int
            default = 100
         "--subsampling","-m"
              help = "Numbers of subgroups, subsampling for stablility selection, pleaase keep default for most users"
              arg_type = Int
              default=2
        "--RNAseqdata","-i"
             help = "Input RNAseq name, required"
             arg_type = String
             required = true
         "--survival_time_status","-t"
                help = "Input time and status file name, required"
                arg_type = String
                required = true
         "--clinic_covariates","-v"
                help = "Input clinic clinic_covariate file name, required"
                arg_type = String
                required = true
         "--pcaname","-p"
                help = "Input PCA score file name, required"
                arg_type = String
                required = true
         "-o"
               help = "Output file name, required"
               arg_type=String
               required =true
    end
    return parse_args(arglist)
end

##read in user input
pargs = userArgs()
addprocs(pargs["cpus"])


function read2array(ifname)
   a=Array{Float32}[]
   geneDic=Dict()
   ifile=open(ifname,"r")
   genes=readline(ifile)
   genes=split(strip(genes),"\t")
   println(genes[1:3])
   [geneDic[i]=genes[i] for i in collect(1:length(genes))]
   for line in eachline(ifile)
        line=split(strip(line),"\t")
        x=[parse(Float32,i) for i in line[1:length(line)]]
        push!(a,x)
   end
  return(geneDic,a)
end



function dicKeys(d)
   dd=Dict()
   for i in 1:length(d)
      dd[string(i)]=d[i]
      for j in i+1:length(d)
         k=string(string(i),"_",string(j))
         dd[k]=string(d[i],"_",d[j])
      end
   end
 return dd
end



@everywhere using RCall
@everywhere R"(library(survival))"
@everywhere using StatsBase
@everywhere using Statistics
@everywhere function coefRcall(a11)
   co_sub=[0 0 0 0]
   rundata=hcat([a11,clinic_covariates]...)
   for i in 1:niteration
      rand_index=sample(collect(1:size(a,1)),size(a,1),replace=false)
      g=Int(ceil(size(a,1)/split_groups))
      split_index=collect(Iterators.partition(rand_index, g)) #
      for j in 1:split_groups
           aj=rundata[split_index[j],:]
           clinicj=clinic[split_index[j],:]
           @rput aj
           @rput clinicj
           co=R"(summary(coxph(Surv(clinicj$time, clinicj$status) ~., data=aj)))"
           rco=rcopy(co)
           jcoef=rco[:coefficients][1,[1,5]]
           jconf=rco[:conf_int][1,[3,4]]
           jcon=[jcoef[1] jcoef[2] jconf[1] jconf[2]]
           co_sub=vcat(co_sub,jcon)
      end
   end
   co_sub=co_sub[2:end,:]
   co_sub=mean(co_sub,dims=(1))
   return co_sub
end


@everywhere function subSet(sub_range)
   sleep(30)
   outfile = open(string(outfilename,"_",string(sub_range[1]),"_",string(sub_range[end]),".txt"), "w")
   sub_dic=Dict()
   for inode in sub_range #run i node
      icoef=coefRcall(a[:,inode])
      if icoef[2] > 0.05
         continue
      end
      println(outfile,inode,"\t",icoef[1],"\t",icoef[2],"\t",icoef[3],"\t",icoef[4])
      sub_dic[inode]=icoef

   end
   close(outfile)
   return sub_dic
end


@everywhere using DelimitedFiles
function write2DelimitedFiles(results)
   aa = Array[]
   push!(aa,["id","coef","pvalue",".95lower","0.95upper"])
   [push!(aa, [genedic[inode],v[1],v[2],v[3],v[4]]) for iDic in results for (inode,v) in iDic]
   outfinalfile = open(string(outfilename,"_final.survival",".txt"), "w")
   writedlm(outfinalfile,aa)
   close(outfinalfile)
end

function writekeys(d)
   ofile = open("allkeys_RNAseq.txt", "w")
   for (k,v) in d
      println(ofile,k,"\t",v)
   end
end


function main()

   println("In running, it might take a while!")
   n=Int(ceil(size(a,2)/cpus))
   x=collect(1:size(a,2))
   run_range=[x[i:min(i+n-1,length(x))] for i in 1:n:length(x)]
   results=pmap((sub_range)->subSet(sub_range),[sub_range for sub_range in run_range])
   write2DelimitedFiles(results)
end



cpus=Int(floor(pargs["cpus"]))
@eval @everywhere cpus=$cpus
outfilename=pargs["o"]
@everywhere outfilename= $outfilename

ifname=pargs["RNAseqdata"]
println(ifname)
genedic,a=read2array(ifname)
writekeys(genedic)

a=hcat(a...)
a=transpose(a)
#a=SharedArray(a[1:end,1:end],)
a = SharedArray{Float32,2}(a[1:end,1:end])
@eval @everywhere a=$a
println("your data, first 2 rows and columns")
println(a[1:2,1:2])

@everywhere using CSV
clinic_cov_name=pargs["clinic_covariates"]
clinic_covariates=CSV.read(clinic_cov_name,delim="\t")
pcaname=pargs["pcaname"]
pca=CSV.read(pcaname,delim="\t")
clinic_covariates=hcat([clinic_covariates,pca]...)
@eval @everywhere clinic_covariates = $clinic_covariates
survival_name=pargs["survival_time_status"]
clinic=CSV.read(survival_name,delim="\t")
@eval @everywhere clinic = $clinic
split_groups =pargs["subsampling"]
@eval @everywhere split_groups = $split_groups
niteration =pargs["iteration"]
@eval @everywhere niteration = $niteration


main()
