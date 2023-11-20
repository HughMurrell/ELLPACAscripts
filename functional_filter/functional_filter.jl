ENV["MPLBACKEND"] = "Agg"
using NextGenSeqUtils, FASTX, BioSequences, StatsBase, SplitApplyCombine, DataFrames, CSV

########### first define some functions for later use

function NextGenSeqUtils.translate_to_aa(s::String)
    s=s[1:3*div(length(s),3)]
    rna = convert(LongRNA{4}, LongDNA{4}(s))
    return string(translate(rna))
end

function reverse_translate_with_gaps(aa_seq::String, nu_seq::String)
    nu_seq = degap(nu_seq)
    if 3*length(degap(aa_seq)) != length(nu_seq)
        println("3 * length of degapped aa_seq = $(3*length(degap(aa_seq)))")
        println(" length of nu_seq = $(length(nu_seq) )")
        println("incompatible lengths")
        return("")
    end
    nu_seq = reverse(collect(nu_seq))
    ret = []
    for ch in collect(aa_seq)
        if ch == '-'
            push!(ret,'-')
            push!(ret,'-')
            push!(ret,'-')
        else
            push!(ret,pop!(nu_seq))
            push!(ret,pop!(nu_seq))
            push!(ret,pop!(nu_seq))
        end
    end
    return(join(ret))
end
        

function my_write_fasta(filename, seqs;
    names=String[], LongSequence = false, append = false, aa = false)
    if !LongSequence
        if aa
            seqs = [BioSequences.LongAA(s) for s in seqs]
        else
            seqs = [BioSequences.LongDNA{4}(s) for s in seqs]
        end
    end
    stream = open(FASTA.Writer, filename, append=append)
    i = 0
    if length(names) != length(seqs)
        names = [string("seq_", i) for i in 1:length(seqs)]
    end
    for (s, n) in zip(seqs, names)
        i += 1
        write(stream, FASTA.Record(n, s))
    end
    close(stream)
end

function read_fasta_with_names_and_descriptions(filename; seqtype=String)
    nams=[]
    seqs=[]
    reader = open(FASTA.Reader, filename)

        for rec in reader
            push!(nams,FASTA.description(rec))
            push!(seqs,FASTA.sequence(String,rec))
        end

    close(reader)
    return nams, seqs
end

"""
    mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)
Julia wrapper for mafft.
"""
function my_mafft(inpath, outpath)
    # cmd = `mafft-fftns --quiet --thread 2 --ep 2 --op 3 --out $outpath $inpath`
    cmd = `mafft --quiet --thread 2 --ep 2 --op 3 --out $outpath $inpath`
    println(cmd)
    run(cmd)
end

function consensus(seqs)
    cons = join([mode([seqs[i][j]
                    for i in 1:length(seqs)])
                        for j in 1:length(seqs[1])])
    return(cons)
end

function agreement(ref,seq)
    if length(ref) == length(seq)
        return( sum(collect(ref).==collect(seq)) / length(ref) )
    else
        return 0.0
    end
end

function get_matching_region_for_consensus(ref, cons)
    # size_start=3, tol_start=0, size_stop=3, tol_stop=0)
    # translate reference
    ref_aa = String(BioSequences.translate(LongDNA{4}(degap(ref))))
    # first generate three reading frames
    cons=degap(cons)
    cons_frames = [ cons[1:3*div(length(cons),3)],
                    cons[2:3*div(length(cons),3)-2],
                    cons[3:3*div(length(cons),3)-1] ]
    # then translate them
    cons_aa_frames = (x->String(BioSequences.translate(LongDNA{4}(degap(x))))).(cons_frames)
    # then choose the frame with the least stop codons
    stop_counts = (x->count('*',x)).(cons_aa_frames)
    # println(stop_counts)
    best_reading_frame = findmin(stop_counts)[2]
    cons_aa = cons_aa_frames[best_reading_frame]
    # find longest coding region in all frames
    start=1
    stop=length(cons_aa)
    max_match_length=0
    for rf in 1:3
        cons_aa = cons_aa_frames[rf]
        m=match(r"M[^\*]*\*",cons_aa)
        while ! isnothing(m)
            if length(m.match)>max_match_length
                best_reading_frame=rf
                start=m.offset
                stop=start-1+length(m.match)
                max_match_length=length(m.match)
            end
            m=match(r"M[^\*]*\*",cons_aa,m.offset+1)
        end
    end
    reading_frame=best_reading_frame
    cons_aa = cons_aa_frames[reading_frame]
    # println(reading_frame)
    # println(start)
    # println(stop)
    # println(cons_aa[start:stop])
    # start_query = ApproximateSearchQuery(LongSequence{AminoAcidAlphabet}(ref_aa[1:size_start]))
    # start = findfirst(start_query, tol_start, LongSequence{AminoAcidAlphabet}(cons_aa))
    # println(cons_aa[start])
    # start=collect(start)[1]
    # cons_aa=cons_aa[start:end]
    # find matching end of ref in translated consensus
    # stop_query = ApproximateSearchQuery(LongSequence{AminoAcidAlphabet}("*"))
    # stop = findfirst(stop_query, 0, LongSequence{AminoAcidAlphabet}(cons_aa))
    # println(stop)
    # println(cons_aa[stop])
    # stop=collect(stop)[end]
    cons_aa_trim = cons_aa[start:stop]
    cons_trim = cons[((start-1)*3+reading_frame):((stop)*3+reading_frame-1)]
    # println(cons_aa_trim)
    return cons_trim, cons_aa_trim
end
                        

function get_matching_region(ref, query; size=10, tol=1, is_cons=false)
    start_query = ApproximateSearchQuery(LongDNA{4}(ref[1:size]))
    start = findfirst(start_query, tol, LongDNA{4}(query))
    # println("start = ",start)
    # stop_query = ApproximateSearchQuery(LongDNA{4}(ref[end-(size*3):end]))
    # stop = findlast(stop_query, tol*4, LongDNA{4}(query))
                                
    stop_query = ApproximateSearchQuery(LongDNA{4}(reverse(ref[end-(size):end])))
    stop = findfirst(stop_query, tol, LongDNA{4}(reverse(query)))
                                
    # println("stop = ",stop)
    if ( isnothing(start) | isnothing(stop) )
        # println("warning --- no matching regions, discarding ... ")
        return(query,translate_to_aa(query),false,"no_matching_region")
    end
    # cr = query[collect(start)[1]:collect(stop)[end]]
    cr = query[collect(start)[1]:length(query)-collect(stop)[1]+1]
    cr_aa = translate_to_aa(cr)
    # now do triplet align of trimmed query to ref and resolve
    # ali = triplet_nw_align(ref, cr, edge_reduction = 0.99, boundary_mult = 2)
    ali = triplet_kmer_seeded_align(ref, cr,
       wordlength = 30, skip = 9, boundary_mult = 2,
       alignedcodons = true, debug=true)
    # println("ali_align = ", ali[2])
    ali = resolve_alignments(ali[1], ali[2], mode = 2)
    cr = ali[2]
    # println("alignment length = ",length(cr))
    # println("cr = ",cr)
    cr = cr[1:(3*div(length(cr),3))]
    # println(cr)
    cr_aa = String(BioSequences.translate(LongDNA{4}(cr)))
    # println("cr_aa = ",cr_aa)
    stops = findall((x->x=='*').(collect(cr_aa)))
    if length(stops) > 0 && stops[1] != length(cr_aa)
        if is_cons
            # cr = cr[1:stops[1]*3]
            # cr_aa = String(BioSequences.translate(LongDNA{4}(cr)))
            # println("warning --- consensus has an internal stops at positions $(stops[1])...")
            return(cr, cr_aa,false,"cons_internal_stop")
        else
            cr = cr[1:stops[1]*3]
            cr_aa = String(BioSequences.translate(LongDNA{4}(cr)))
            # println("warning --- internal stop codon at $(stops[1]), discarding ... ")
            return(cr,cr_aa,false,"internal_stop")
        end
    end
    frame_errs =  findall((x->x=='X').(collect(cr_aa)))
    if length(frame_errs) > 0
        if is_cons
            # cr = cr[1:stops[1]*3]
            # cr_aa = String(BioSequences.translate(LongDNA{4}(cr)))
            # println("warning --- consensus has frame error...")
            return(cr,cr_aa,false,"cons_frame_error")
        else
            # println("warning --- frame error, discarding ... ")
            return(cr,cr_aa,false,"frame_error")
        end
    end
    if (cr_aa[1] != 'M') | (cr_aa[end] != '*')
        # println("warning --- seq has no START-and-STOP codon, discarding ... ")
        return(cr,cr_aa,false,"no_start_stop")
    end
    return(cr, cr_aa, true, "")
end



function extract_ref_match(ref_file, donor, seq_files ; house_keeping = nothing )
    
    # set the output path names
    faa=mkpath("data_out/functional_aa/")
    out_seq_aa = faa*"/"*donor*"_aa_sequences.fasta"
    faa_ali=mkpath("data_out/functional_aa_ali/")
    out_ali_aa = faa_ali*"/"*donor*"_aa_alignment.fasta"
    fnu=mkpath("data_out/functional_nu/")
    out_seq_nu = fnu*"/"*donor*"_nu_ca_sequences.fasta"
    fnu_ali=mkpath("data_out/functional_nu_ali/")
    out_ali_nu = fnu_ali*"/"*donor*"_nu_ca_alignment.fasta"
    nfaa=mkpath("data_out/non_functional_aa/")
    nf_out_seq_aa = nfaa*"/"*donor*"_aa_sequences_non_functional_rejects.fasta"
    nfnu=mkpath("data_out/non_functional_nu/")
    nf_out_seq_nu = nfnu*"/"*donor*"_nu_ca_sequences_non_functional_rejects.fasta"
    
    # first read and translate reference file
    println("processing ref and sequence files using NextGenSeqUtils...")
    println("reading reference...")
    nr, sr = read_fasta(ref_file)
    sr_aa = String(BioSequences.translate(LongDNA{4}(degap(sr[1]))))
    if (sr_aa[1] != 'M') | (sr_aa[end] != '*')
        error("bad reference")
        exit()
    end
    my_write_fasta(out_seq_aa,[sr_aa],names=[nr[1]],aa=true,append=false)
    my_write_fasta(nf_out_seq_aa,[sr_aa],names=[nr[1]],aa=true,append=false)
    my_write_fasta(out_seq_nu,[degap(sr[1])],names=[nr[1]],aa=false,append=false)
    my_write_fasta(nf_out_seq_nu,[degap(sr[1])],names=[nr[1]],aa=false,append=false)
    println("reference read, translated and written....")
    
    # get consensus of first sequence file and prepare a pannel sequence
    file = seq_files[1]
    println("processing first sequence file with name $(basename(file)) ...")
    ns, ss = read_fasta_with_names_and_descriptions(file)
    cons = uppercase(degap(consensus(ss)))
    println("length of reference = $(length(sr[1]))")
    println("length of consensus = $(length(cons))")
    #cons, cons_aa = get_matching_region(sr[1], cons, size=9, tol=1, is_cons=true)
    cons, cons_aa = get_matching_region_for_consensus(sr[1], cons)
    # println(cons)
    # println(cons_aa)
    stops = findall((x->x=='*').(collect(cons_aa)))
    if length(stops) == 0 || stops[1] != length(cons_aa)
        # try a direct search for coding region
        cons = uppercase(degap(consensus(ss)))
        cons,cons_aa=get_matching_region_for_consensus(sr[1], cons)
    end
    println("length of matching region = $(length(cons))")
    if length(cons)==0
        println("Consensus has no matching regions...")
        return()
    end
    println("consensus generated and coding region extracted...")
    my_write_fasta(out_seq_aa,[cons_aa],
        names=["consensus_of_first_sample"], aa=true,append=true)
    # also write the nucs for later use
    my_write_fasta(out_seq_nu,[degap(cons)],
        names=["consensus_of_first_sample"], aa=false, append=true)

    # now read all sequences and extract reading frames by aligning to cons
    # println("now extracting coding regions from seqs by comparing to consensus...")
    # println()
    for k in 1:length(seq_files)
        file = seq_files[k]
        pool=match(r"pool.",file).match[5]
        ns, ss = read_fasta_with_names_and_descriptions(file)
        for j in 1:length(ns)
            splits=split(ns[j]," ")
            nam1=splits[1]*"_pool$(pool)"
            desc=""
            if length(splits) > 1
                desc=join(splits[2:end]," ")
            end
            ns[j]=join([nam1,desc]," ")
        end
        println("processing $(length(ss)) sequences from file $(basename(file))  ")
        count_good = 0
        count_bad = 0
        count_nmr = 0
        count_isc = 0
        count_rfe = 0
        count_nss = 0
        count_length = 0
        for i in 1:length(ss)
            ts = uppercase(degap(ss[i]))
            ts, ts_aa, ok, msg = get_matching_region(cons, ts, size=15, tol=2)
            # ts, ts_aa = get_matching_region_for_consensus(cons, ts)
            # println("length of ts_aa = $(length(ts_aa))")
            if ok
                count_length += length(ts_aa)
                my_write_fasta(out_seq_aa,[ts_aa],
                    names=[ns[i]], aa=true, append=true)
                my_write_fasta(out_seq_nu,[degap(ts)],
                    names=[ns[i]], aa=false, append=true)
                count_good += 1
            else
                my_write_fasta(nf_out_seq_aa,[ts_aa],
                    names=[ns[i]*" "*msg], aa=true, append=true)
                my_write_fasta(nf_out_seq_nu,[degap(ts)],
                    names=[ns[i]*" "*msg], aa=false, append=true)
                count_bad += 1
                msg == "no_matching_region" ? count_nmr += 1 : nothing
                msg == "internal_stop" ? count_isc += 1 : nothing
                msg == "frame_error" ? count_rfe += 1 : nothing
                msg == "no_start_stop" ? count_nss += 1 : nothing
            end
        end
        print("processed $(count_good) good seqs and $(count_bad) bad seqs, ")
        println("success rate = $(count_good / (count_good + count_bad))")
        if ! isnothing(house_keeping)
            visit=basename(file)[8:11]
            mean_aa_length = count_good > 0 ? Int(round(count_length / count_good)) : 0
            push!(house_keeping,[donor,visit,pool,count_good,count_bad,
                        mean_aa_length,count_nmr,count_isc,count_rfe,count_nss])
        end
    end
    println()
    println("using mafft to align aa sequences (please be patient)....")
    my_mafft(out_seq_aa, out_ali_aa)
    
    # println("Doing min cover checks ...")
    # nr, sr = read_fasta_with_names_and_descriptions(out_ali_aa)
    # ok_inds = ( (x->agreement(sr[1],x)).(sr) ) .>= min_cover
    # if ! ok_inds[2]
    #     println("ERROR: consensus outside minimum cover $(min_cover)")
    # end
    # println("discarding $(length(sr) - sum(ok_inds)) seqs of $(length(sr)) total seqs due to min cover < $(min_cover)...")
    # println("success rate = $(sum(ok_inds) / length(sr))")
    # sr_keeps = sr[ok_inds]
    # nr_keeps = nr[ok_inds]
    # my_write_fasta(out_ali_aa, sr_keeps, names=nr_keeps, aa=true, append=false)

    # sr_discards = sr[(!).(ok_inds)]
    # nr_discards = nr[(!).(ok_inds)]
    # nr_discards = (s->s*" failed_min_cover_check").(nr_discards)
    # my_write_fasta(nf_out_seq_aa, sr_discards, names=nr_discards, aa=true, append=true)
    
    println("generating codon aware nucleotide alignment ...")
    nr, sr = read_fasta_with_names_and_descriptions(out_ali_aa)
    nu_n, nu_s = read_fasta_with_names_and_descriptions(out_seq_nu)
    # nu_n = nu_n[ok_inds]
    # nu_s = nu_s[ok_inds]
    rm(out_ali_nu,force=true)
    for i in 1:length(nu_s)
        re_gapped_ns = reverse_translate_with_gaps(sr[i]*"*",nu_s[i])
        my_write_fasta(out_ali_nu,[re_gapped_ns], names=[nu_n[i]], aa=false, append=true)
    end
    
    println("sorting rejects on reason .....")
    nu_n, nu_s = read_fasta_with_names_and_descriptions(nf_out_seq_nu)
    my_write_fasta(nf_out_seq_nu,[nu_s[1]], names=[nu_n[1]], aa=false, append=false)
    if length(nu_n) > 1
        nu_n=nu_n[2:end]
        nu_s=nu_s[2:end]
        reasons=(s->split(s," ")[4]).(nu_n)
        p=sortperm(reasons)
        nu_n=nu_n[p]
        nu_s=nu_s[p]
        my_write_fasta(nf_out_seq_nu,nu_s, names=nu_n, aa=false, append=true)
    end
    aa_n, aa_s = read_fasta_with_names_and_descriptions(nf_out_seq_aa)
    my_write_fasta(nf_out_seq_aa,[aa_s[1]], names=[aa_n[1]], aa=true, append=false)
    if length(aa_n) > 1
        aa_n=aa_n[2:end]
        aa_s=aa_s[2:end]
        reasons=(s->split(s," ")[4]).(aa_n)
        p=sortperm(reasons)
        aa_n=aa_n[p]
        aa_s=aa_s[p]
        my_write_fasta(nf_out_seq_aa,aa_s, names=aa_n, aa=true, append=true)
    end
    println("$(donor) completed ...")
end

############## end of function definitions


println("using Julia version: $(VERSION)")
t1 = time()
num_of_pools=4
pool_dirs=["../../../porpidpostproc_ellp/postproc/pool$(i)_" for i in 1:num_of_pools]
for i in 1:num_of_pools
    if i<3
        pool_dirs[i] *= "nicd"
    elseif i<4
        pool_dirs[i] *= "cell123"
    elseif i<5
        pool_dirs[i] *= "earlham_mark2"
    end
end
# println(pool_dirs)

pool_sample_dict=Dict()
for i in 1:num_of_pools
    samples = readdir(pool_dirs[i])
    samples = samples[(s -> ~ contains(s,'.')).(samples)]
    println("$(pool_dirs[i]) ======> $(length(samples)) samples")
    pool_sample_dict[i]=samples
end

donors=[]
for k in keys(pool_sample_dict)
    # println("$(k) -> $(length(pool_sample_dict[k])) samples")
    pd=(s -> s[1:6]).(pool_sample_dict[k])
    global donors = vcat(donors,pd)
end

donors=sort(union(donors))
# println(donors)
println("Number of donors = $(length(donors))")

fr="../../../porpidpostproc_ellp/panels/hxb2-env.fasta"
println(" doing codon aware alignment against reference in $(fr) ...")

hk=DataFrame(donor=[],visit=[],pool=[],functional=[],non_functional=[],
                        mean_functional_aa_length=[],
                        no_matching_region=[],internal_stop_codon=[],
                        reading_frame_error=[],no_stop_and_start=[])
for donor in donors
    # func_path = mkpath("functional/"*donor)*"/"
    # non_func_path = mkpath("non_functional/"*donor)*"/"
    fs = []
    for pool in 1:num_of_pools
        for sample in pool_sample_dict[pool]
            if sample[1:6]==donor
                push!(fs,pool_dirs[pool]*"/"*sample*"/"*sample*".fasta")
            end
        end
    end
    si=(x->split(x,"/")[7]).(fs)
    fs=fs[sortperm(si)]
    println()
    println("------- $(donor) -------")
    extract_ref_match(fr, donor, fs, house_keeping=hk)
end

println("writing housekeeping report ...")
mkpath("reports/")
CSV.write("reports/housekeeping.csv",hk)

t2 = time()
println("Functional sequences filtered in $(t2 - t1) seconds.")
