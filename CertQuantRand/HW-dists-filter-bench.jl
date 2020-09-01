### Luís Brandão (NIST/Strativia), 2020-September-01
### Programming code written in Julia language v1.3

__precompile__()

module mod_HWdistFilters

### DESCRIPTION

# Module to benchmark the filtering of low hamming-weight (HW) distances (dists)

#= Consider a set with k strings (each with n bits, for some n<=64), possibly
   with repetitions. The goal is to apply a filter that outputs a subset of
   strings such that any pair of strings has HW-dist higher than a threshold t.
   For example, in the context of certifiable quantum-randomness experiments,
   one may be interested in parameters like (k,n,t)=(10^6,53,10).

   As currently imlemented, for practical parameters of interest, the bottleneck
   is the check (with complexity quadratic in k) of the HW-dist across all
   possible pairs of strings. Once the connected edges are identified, the
   part of deciding which ones to exclude (i.e., applying a filter) is much
   faster (both in the Linear and in the Greedy approahes).

   Two filters are experimented: Linear and greedy.

   === Linear
   Iterate over all available strings, using a random ordering (e.g., the
   ordering by which strings were initially sampled). For each selected string,
   exclude (i.e., mark as unavailable) all neighbors (i.e., with HW-dist <=t).
   When the iteration is finished, output the subset of selected strings.

   === Greedy
   Iterate over all available strings with positive connectivity (i.e., all
   strings that have neighbors), using a descending ordering of connectivity
   degree. Remove each selected string, and correspondingly update the
   connectivity of its neighbors. Once all remaining strings have connectivity
   degree 0, output he set of those strings.

=#

    export fGenVecStrings
    export fCreateEncodedEdges
    export fParseArrayToEdges
    export fGenDictionariesFromEdges
    export fRemoveNodesFilterGreedy
    export fRemoveNodesFilterLinear
    export fInternalBenchHWdistFilters
    export fBenchHWdistFilters

    using Dates  # function now()

    """
        fGenVecStrings(n::Int64,k::Int64)::Vector{UInt64}
    ## INPUT
        - `n`::Int64: number of significant bits per string (max allowed is 64)
        - `k`::Int64: number of strings to generate
    """
    function fGenVecStrings(n::Int64,k::Int64)::Vector{UInt64}
        if n>64; error("Code only works with n up to 64 qubits") end
        listR=rand(UInt64,k)
        listR .>>>= (64-n)
    end


    """
    fCreateEncodedEdges(listNodes::Vector{Int64},thresholdHW::Int64)
    ##INPUT
        - `listNodes`::Vector{Int64}: list of UInt64s, each representing a
           string, where only the n LSBs are relevant\n
        - `thresholdHW`:Int64: threshold of HW distance (must be <= n)\n
    ##OUTPUT
        - encodedListEdges::Vector{Int64}: an intermediate representation of a
        list of edges, encoded as a simple vector of Int64 indices. A value 0
        indicates incrementing the first coordinate of an edge. An edge is
         created between two nodes iff their HW distance is <= thresholdHW.
        Later, an edge will be a pair of indices
    """
    function fCreateEncodedEdges(listNodes,thresholdHW::Int64)::Vector{Int64}
        k::Int64=length(listNodes)
        len::Int64=2*k
        z::Int64=0 #z: position in encodedListEdges
        encodedListEdges::Vector{Int64}=zeros(Int64,len)
        L2::Int64=len-k
        for i=1:k;
            s1=listNodes[i]; z+=1; #encodedListEdges[z]=0 # It's already 0
            # Could avoid the next `if' if predicting `len' with a good margin
            # len+=2k; append!(encodedListEdges,zeros(Int64,2k)); end
            if L2 < z; len+=2k;   # same as `if len-z < k', but no subtraction
                append!(encodedListEdges,zeros(Int64,2k)); L2=len-k; end
            for j::Int64=i+1:k; s2=listNodes[j]
                dist=count_ones(xor(s1,s2))
                if dist<=thresholdHW; z+=1; encodedListEdges[z]=j end
            end
        end
        encodedListEdges[1:z]
    end


    """
    fParseArrayToEdges(encodedListEdges::Vector{Int64})
    ##INPUT
        - `encodedListEdges`::Vector{Int64}: as outputted by
          fCreateEncodedEdges(listNodes,thresholdHW)
    ##OUTPUT
        - listEdges::Vector{Tuple{Int64,Int64}: a list of pairs, where
          each pair is an edge connecting two indices
    """
    function fParseArrayToEdges(encodedListEdges::Vector{Int64})::Vector{Tuple{Int64,Int64}}
        i1=0
        listEdges=Vector{Tuple{Int64,Int64}}()
        for i=1:length(encodedListEdges)
            i2=encodedListEdges[i]
            if i2==0; i1+=1
            else push!(listEdges,(i1,i2))
            end
        end
        listEdges
    end


    """fGenDictionariesFromEdges(listEdges)
    ##INPUT
        - k::Int64: number of possible strings (i.e., of possible indices)\n
        - listEdges::Vector{Tuple{Int64,Int64}}: each edge is (i1,i2),
          where 1<=i1,i2<=k
    ##OUTPUT
        A quartet (maxDegree,dCount,dConnect,dDegreeSet) of dictionaries:\n
        - maxDegree::Int64: the maximum degree of a node\n
        - dCount::Dict{Int64,Int64}(): dCount[i1]=n1 means that there are n1
          strings with HW distance <= t from string i1
        - dConnect::Dict{Int64,Vector{Int64}}(): dCount[i1]=[n1,n2,...] means
          that [n1,n2,...] are the strings at a distance <= t from i1
        - dDegreeSet::Dict{Int64,Set{Int64}}(): dCount[d1]=Set(i1,i2,...) means
          that each of the strings i1,i2... have connectivity degree d1
    """
    function fGenDictionariesFromEdges(k::Int64,
                listEdges::Vector{Tuple{Int64,Int64}})::Tuple{Int64,
                Dict{Int64,Int64},Dict{Int64,Set{Int64}},Dict{Int64,Set{Int64}}}
        # dCount: dictionary: for each string, counts how many neighbors.
        dCount=Dict{Int64,Int64}();
        for i=1:k; dCount[i]=0; end
        # dConnect: dictionary: for each string, stores the set of neighbors.
        dConnect=Dict{Int64,Set{Int64}}();
        for i=1:k; dConnect[i]=Set{Int64}(); end
        # dDegreeSet: for each degree d, maps the strings with connectivity d.
        dDegreeSet=Dict{Int64,Set{Int64}}();
        ###
        for (i1,i2) in listEdges
            dCount[i1]+=1; push!(dConnect[i1],i2);
            dCount[i2]+=1; push!(dConnect[i2],i1)
        end
        maxdeg=0
        for i=1:k
            deg=dCount[i]
            if deg>maxdeg; maxdeg=deg end
            if !haskey(dDegreeSet,deg)
                dDegreeSet[deg]=Set{Int64}([i])
            else push!(dDegreeSet[deg],i)
            end
        end
        # Initializes as empty vector the unitialized cases of degree < maxdeg
        for i=0:maxdeg-1;
            if !haskey(dDegreeSet,i); dDegreeSet[i]=Set{Int64}() end
        end
        (maxdeg,dCount,dConnect,dDegreeSet)
    end
    1+1


    """
        fRemoveNodesFilterGreedy(maxdeg,dCount1,dConnect1d,DegreeSet1)
    ##INPUT:
        maxdeg,dCount,dConnect,dDegreeSet (output by fGenDictionariesFromEdges)
    ##OUTPUT
        list of string indices, corresponding to the list of strings
        (e.g., as previously output by fGenVecStrings)
    """
    function fRemoveNodesFilterGreedy(maxdeg::Int64,
                            dCount1::Dict{Int64,Int64},
                            dConnect1::Dict{Int64,Set{Int64}},
                            dDegreeSet1::Dict{Int64,Set{Int64}})::Vector{Int64}
        ### Need to deepcopy the dictionaries because they will be changed
        dCount=deepcopy(dCount1)
        dConnect=deepcopy(dConnect1)
        dDegreeSet=deepcopy(dDegreeSet1)
        if maxdeg<1; return dDegreeSet[0] end
        for deg=maxdeg:-1:1
            while !isempty(dDegreeSet[deg])
                i=pop!(dDegreeSet[deg])
                for j in dConnect[i];
                    if !in(i,dConnect[j]);
                        error("i=$i was expected inside dConnect[$j]")
                    end
                    #Updates due to node i having been removed
                    delete!(dConnect[j],i)
                    deg2=dCount[j]
                    delete!(dDegreeSet[deg2],j)
                    dCount[j]-=1
                    push!(dDegreeSet[deg2-1],j)
                end
            end
        end
        ### Sanity check
        for deg=1:maxdeg; if !isempty(dDegreeSet[deg]);
            error("Non empty dDegreeSet[$deg]") end end
        return [i for i in dDegreeSet[0]]
    end


    """
        fRemoveNodesFilterLinear(k,dConnectd)
    ##INPUT:
        k: original number of strings (indexed between 1 and k)
        dConnect: as output by fGenDictionariesFromEdges
    ##OUTPUT
        list of string indices, corresponding to the list of strings
        (e.g., as previously output by fGenVecStrings)
    """
    function fRemoveNodesFilterLinear(k::Int64,
                                dConnect::Dict{Int64,Set{Int64}})::Vector{Int64}
        # Boolean marking whether a string is on or not
        listBoolsStrings::Vector{Bool}=trues(k); # BitVector::=repeat([true],k)
        listFinalNodes=Vector{Int64}()
        for i=1:k; if listBoolsStrings[i]
            for j in dConnect[i]; listBoolsStrings[j]=false end
            push!(listFinalNodes,i)
        end end
        listFinalNodes #indices
    end



    function fInternalBenchHWdistFilters(n::Int64,listNodes::Vector{UInt64},
                thresh::Int64)::Tuple{Int64,Int64}
        print("  Encoded edges: ");
        @time encodedListEdges::Vector{Int64}=fCreateEncodedEdges(listNodes,thresh);
        print("  Parse selected edges: ");
        @time listEdges=fParseArrayToEdges(encodedListEdges)
        print("  Connectivity structure: ");
        k=length(listNodes)
        @time (maxdeg,d1,d2,d3)=fGenDictionariesFromEdges(k,listEdges)
        (dCount,dConnect,dDegreeSet)=(d1,d2,d3)
        print("  Apply filter Linear: ");
        @time L1=fRemoveNodesFilterLinear(k,dConnect)
        print("  Apply filter Greedy: ");
        @time L2=fRemoveNodesFilterGreedy(maxdeg,dCount,dConnect,dDegreeSet)
        (k1,k2)=map(length,(L1,L2))
    end


    """
        fBenchHWdistFilters(n,k,thresh)
    ##INPUT:
        - n:
        - k: original number of strings (indexed between 1 and k)
        - thresh:
    ##OUTPUT
        (k1::Int64,k2::Int64,tim::Float64)\n
        - k1 and k2: number of strings (and overall computation time) remaining
        after the Linear and the Greedy filtering, respectively.
        - tim: time (rounded to 4 decimals) for the overall computation
    """
    function fBenchHWdistFilters(n::Int64,k::Int64,
                                    thresh::Int64)::Tuple{Int64,Int64,Float64}
        fName="fBenchHWdistFilters"
        tim1=time(); println("\nBEGIN $fName($n,$k,$thresh): ",now());
        println("  Number of bits per string: n=$n");
        println("  Number of original strings: k=$k");
        println("  Threshold HW-dist of connectivity: t=$thresh");
        print("  Generate random strings: ");
        @time listNodes=fGenVecStrings(n,k)
        (k1,k2)=fInternalBenchHWdistFilters(n,listNodes,thresh)
        tim2=time();
        print("  Remaining strings: ");
        p1=round(k1/k*100,digits=6)
        p2=round(k2/k*100,digits=6)
        println("$p1% ($k1) (Linear); $p2% ($k2) (Greedy).")
        timeDur=round(tim2-tim1,digits=4)
        println("END $fName: ",now(),". Duration: $timeDur seconds");
        (k1,k2,timeDur)
    end

end  ### END OF MODULE mod_HWdistFilters

################################################################################
################################################################################


using .mod_HWdistFilters

### Quick benchmarking ... the first run may be slower
(n,k,thresh)=(53,10^5,10); listNodes=fGenVecStrings(n,k); @time fCreateEncodedEdges(listNodes,thresh);

### Full benchmarking
(n,k,thresh)=(53,10^4,10); fBenchHWdistFilters(n,k,thresh)
(n,k,thresh)=(53,3*10^4,10); fBenchHWdistFilters(n,k,thresh)
(n,k,thresh)=(53,10^5,10); fBenchHWdistFilters(n,k,thresh)
(n,k,thresh)=(53,3*10^5,10); fBenchHWdistFilters(n,k,thresh)
(n,k,thresh)=(53,10^6,10); fBenchHWdistFilters(n,k,thresh)
