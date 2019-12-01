module JuliaChildsGame

greet() = print("Hello World!")
help() = println("Get_Art_Seq(t = # sequences, motif string, l = length of seq, # of mutations to perform to motif)\n", "callGame(sequences in string format, length of motif, number of sequences, number of generations, mutation rate (0.01, 0.001 ...), PWN = find secondary motifs True/False)\n")

function Get_Art_Seq(t, motif, l, mutations)
    #41% of l = number of G's or C's to have
    tmp = (41*l)/100
    locations = []
    #mutations = rand(1:length(motif))
    #println("num mutations", mutations)
    #println("motif ", motif)
    motifs = []
    ATGC = ['A', 'T', 'G', 'C']
    if mutations != 0
        for i in 1:t
            temp = deepcopy(motif)
            #println(temp)
            for count in 1:mutations
                position = rand(1:length(motif))
                #random = rand([0,1])
                #println("random is ", random)
                #if random == 0
                    #println("changing pos")
                    #random_letter = rand(1:4)
                    #holder = string(ATGC[random_letter])
                    #println(typeof(holder))
                    #println(holder)
                    #println(typeof(temp))
                tmpmotif = collect(temp)
                ranBP = rand("ACGT")
                while ranBP == tmpmotif[position]
                    ranBP = rand("ACGT")
                end
                tmpmotif[position] = ranBP
                temp = join(tmpmotif)
                #end
            end
            push!(motifs, temp)
        end
    else
        for i in 1:t
            push!(motifs, motif)
        end
    end
    #println(motifs)
    #print(tmp, "\n")
    #round up
    gcContentNeeded = ceil(Int, tmp)
    #print(gcContentNeeded)
    #countCG = number of C's and G's in motif
    countCG_start = 0
    for i in 1:length(motif)
        if (motif[i] == 'C') || (motif[i] == 'G')
            countCG_start = countCG_start + 1
        end
    end
    #print(gcContentNeeded)
    DNA_SEQ = []
    #for the number of seq needed
    for i in 1:t
        #random i position to start motif as long as i doesn't = too close to end of length
        random_motif_pos = rand(1:(l-(length(motifs[i]))))
        push!(locations, random_motif_pos)
        #for length of seq
        temp = ""
        countCG = countCG_start
        #print(countCG, "\n")
        #print(length(motif), "\n")
        for x in 1:(l-length(motifs[i])+1)
            #if at posiiton i, then add in motif
            if (x == 1) && (x != random_motif_pos)
                #println("case1")
                random_letter = rand(1:4)
                temp = string(temp, "", ATGC[random_letter])
                if (temp[x] == 'G') || (temp[x] == 'C')
                    countCG = countCG + 1
                end
            elseif (x == random_motif_pos) && (x == 1)
                #println("case2")
                temp = deepcopy(motifs[i])
                #update length counter so don't make sequence motif - 1 too big
                #print(x, "\n")
                x = x + (length(motifs[i]))
                #print(x, "\n")
            elseif (x == random_motif_pos) && (x > 1)
                #println("case3")
                temp = string(temp, "", motifs[i])
                #print(x, "\n")
                x = x + (length(motifs[i]))
                #print(x, "\n")
            #if countCG >= 41% of length 
            elseif (countCG == gcContentNeeded)
                #println("case4")
                #only pick from A or T
                random_letter = rand(1:2)
                temp = string(temp, "", ATGC[random_letter])
            #if (l - countCG) >= (100%-41% of length)
            elseif ((l - countCG) == (1-gcContentNeeded))
                #println("case5")
                #only pick from G or C
                random_letter = rand(3:4)
                temp = string(temp, "", ATGC[random_letter])
                if (temp[x] == 'G') || (temp[x] == 'C')
                    countCG = countCG + 1
                end
            #else
            else
                #println("case6")
                #randomly pick from ACGorT
                random_letter = rand(1:4)
                temp = string(temp, "", ATGC[random_letter])
                if (temp[x] == 'G') || (temp[x] == 'C')
                    countCG = countCG + 1
                end
            end
            #println(temp)
        end
        #print(countCG, "\n")
        if i == 1
            DNA_SEQ = temp
        else
            DNA_SEQ = hcat(DNA_SEQ, temp)
        end
    end
    #println(DNA_SEQ)
    #println(motifs)
    #println(locations)
    return(DNA_SEQ)
end

# count matrix
function bpCount(motifs, t, k)
    Acount = zeros(1,k)
    Ccount = zeros(1,k)
    Gcount = zeros(1,k)
    Tcount = zeros(1,k)
    for i in 1:k
        for j in 1:t
            if motifs[j][i] == 'A'
                Acount[i] = Acount[i] + 1
            elseif motifs[j][i] == 'C'
                Ccount[i] = Ccount[i] + 1
            elseif motifs[j][i] == 'G'
                Gcount[i] = Gcount[i] + 1
            else
                Tcount[i] = Tcount[i] + 1
            end
        end
        Acount[i] = Acount[i] + 1
        Ccount[i] = Ccount[i] + 1
        Gcount[i] = Gcount[i] + 1
        Tcount[i] = Tcount[i] + 1
    end
    return(Acount, Ccount, Gcount, Tcount)
end

# motif freq matrix
function motFreq(Acount, Ccount, Gcount, Tcount, t, k)
    Amotfreq = zeros(1,k)
    Cmotfreq = zeros(1,k)
    Gmotfreq = zeros(1,k)
    Tmotfreq = zeros(1,k)
    for i in 1:k
        totalcount = Acount[i] + Ccount[i] + Gcount[i] + Tcount[i]
        Amotfreq[i] = Acount[i]/totalcount
        Cmotfreq[i] = Ccount[i]/totalcount
        Gmotfreq[i] = Gcount[i]/totalcount
        Tmotfreq[i] = Tcount[i]/totalcount
    end
    return(Amotfreq, Cmotfreq, Gmotfreq, Tmotfreq)
end

# background counts
function backgroundCounts(backgrounds, t, k)
    Atotcount = 0
    Ctotcount = 0
    Gtotcount = 0
    Ttotcount = 0
    for i in 1:t
        for j in 1:length(backgrounds[i])
            if backgrounds[i][j] == 'A'
                Atotcount = Atotcount + 1
            elseif backgrounds[i][j] == 'C'
                Ctotcount = Ctotcount + 1
            elseif backgrounds[i][j] == 'G'
                Gtotcount = Gtotcount + 1
            else
                Ttotcount = Ttotcount + 1
            end
        end
    end
    # background frequencies
    totaltotalcount = Ttotcount + Atotcount + Ctotcount + Gtotcount
    Atotfreq = Atotcount/totaltotalcount
    Ctotfreq = Ctotcount/totaltotalcount
    Gtotfreq = Gtotcount/totaltotalcount
    Ttotfreq = Ttotcount/totaltotalcount
    return(Atotfreq, Ctotfreq, Gtotfreq, Ttotfreq)
end

# calculating the score function
function product(motifs, backgrounds, t, k)
    Acount, Ccount, Gcount, Tcount = bpCount(motifs, t, k)
    Atotfreq, Ctotfreq, Gtotfreq, Ttotfreq = backgroundCounts(backgrounds, t, k)
    Amotfreq, Cmotfreq, Gmotfreq, Tmotfreq = motFreq(Acount, Ccount, Gcount, Tcount, t, k)
    productthing = 1
    for i in 1:k
        for j in 1:4
            if j == 1
                thetajk = Amotfreq[i]
                thetaok = Atotfreq
            elseif j == 2
                thetajk = Cmotfreq[i]
                thetaok = Ctotfreq
            elseif j == 3
                thetajk = Gmotfreq[i]
                thetaok = Gtotfreq
            else
                thetajk = Tmotfreq[i]
                thetaok = Ttotfreq
            end 
            logthing = log(thetajk/thetaok)

            tempproduct = thetajk*logthing

            productthing = productthing*tempproduct

        end
    end

    totl = 0
    for i in 1:t
        totl = totl + length(backgrounds[i]) + 1
    end

    phat = t/totl

    secondlogthing = log(phat/(1-phat))

    scoremotif = (t*(secondlogthing - 1 + productthing))

    productthing
    return(productthing)
end 

function mutation(motifs, t, r, k, dna)
    randomnumber2 = rand()
    check = false
    newMotifSeq = 0
    tempmotif = 0
    tempback = 0
    if randomnumber2 <= r
        #print(motifs, "\n")
        randomseq = rand(1:t)
        randpos = rand(1:length(dna[randomseq])-k+1)
        tempmotif = SubString(dna[randomseq],randpos,randpos+k-1)
        #motifs[randomseq] = tempmotif
        #backgrounds[randomseq] = tempback
        if randpos == 1
        
        prefix = ""
        suffix = SubString(dna[randomseq],randpos+k,length(dna[randomseq]))
        
        elseif randpos == length(dna[randomseq])-k+1
        prefix = SubString(dna[randomseq],1,randpos-1)
        suffix = ""
        
        else
        prefix = SubString(dna[randomseq],1,randpos-1)
        suffix = SubString(dna[randomseq],randpos+k,length(dna[randomseq]))
        end
        newMotifSeq = randomseq - 1
        check = true
    end
    return (check, newMotifSeq, tempmotif, tempback)
end

function adjust(motifs_A, backgrounds_A, t, k, dna)
    global new_back = ""
    global output = ""
    max_score = 0
    for m in 1:100
        #println(m)
        random = rand(1:t)
        #println(random)
        DNA_Sequence = deepcopy(motifs_A[random])
        motifs_A_notI = deepcopy(motifs_A[1:end .!= random])
        backgrounds_A_notI = deepcopy(backgrounds_A[1:end .!= random])
        Acount, Ccount, Gcount, Tcount = bpCount(motifs_A_notI, t-1, k)
        Amotfreq, Cmotfreq, Gmotfreq, Tmotfreq = motFreq(Acount, Ccount, Gcount, Tcount, t-1, k)
        #println(Amotfreq)
        #println(Cmotfreq)
        #println(Gmotfreq)
        #println(Tmotfreq)
        global output = ""
        prob_start = 1
        for j in 1:length(DNA_Sequence)
                if DNA_Sequence[j] == 'A'
                    prob_start = prob_start*Amotfreq[j]
                elseif DNA_Sequence[j] == 'C'
                    prob_start = prob_start*Cmotfreq[j]
                elseif DNA_Sequence[j] == 'G'
                    prob_start = prob_start*Gmotfreq[j]
                else
                    prob_start = prob_start*Tmotfreq[j]
                end
        end
        #println(prob_start)
        if m == 1
            max_score = prob_start
        end
        score1 = product(motifs_A_notI, backgrounds_A_notI, t-1, k)
        temp_backgrounds = ""
        subsetvec = []
        probvec = []

        #new_motif_location = A vector that is number of DNA sequences
        i = 1
        while i < length(dna[random])-k+1
            #println(i) #SubString(dna[i],i,i+k-1)
            subsetDNA = SubString(dna[random],i,i+k-1)
            #println(subsetDNA)
            prob = 1
            for j in 1:length(subsetDNA)
                if subsetDNA[j] == 'A'
                    prob = prob*Amotfreq[j]
                elseif subsetDNA[j] == 'C'
                    prob = prob*Cmotfreq[j]
                elseif subsetDNA[j] == 'G'
                    prob = prob*Gmotfreq[j]
                else
                    prob = prob*Tmotfreq[j]
                end
            end
                if i == 1
                    prefix = ""
                    suffix = SubString(dna[random],i+k,length(dna[random]))
                elseif i == length(dna[random])-k+1
                    prefix = SubString(dna[random],1,i-1)
                    suffix = ""
                else
                    prefix = SubString(dna[random],1,i-1)
                    suffix = SubString(dna[random],i+k,length(dna[random]))
                end
                temp_backgrounds = string(prefix,"",suffix)
            push!(subsetvec, subsetDNA)
            push!(probvec, prob)
            #println(prob)
            #println(maximum(probvec))
            if i == 1 && prob > max_score
                global output = subsetDNA
                global new_back = temp_backgrounds
                max_score = prob
                #println(output)
            else
                if prob == maximum(probvec) && prob > max_score
                    global output = subsetDNA
                    global new_back = temp_backgrounds
                    max_score = prob
                    #println(output)
                    #new_motif_location[random_number] = i                
                end
            end

            i = i + 1
        end
        #println(random)
        #println("-")
        if output != ""
            motifs_A[random] = output
            backgrounds_A[random] = new_back
        end
        #println(motifs_A)
            


        ##un this a million times over

    end
    #println(motifs_A)
    return(motifs_A, backgrounds_A)
end

function shift(motifs_B, backgrounds_B, t, k, dna)
    #shift all motifs in set over by 1
    temp = motifs_B
    temp_backgrounds = backgrounds_B
    for j in 1:5
        score1 = product(motifs_B, backgrounds_B, t, k)
        for i in 1:t
            #return 1:2
            temp = motifs_B
            temp_backgrounds = backgrounds_B
            #println(findfirst(motifs_B[i], dna[i])[1]+ 1)
            position = findfirst(motifs_B[i], dna[i])[1] + 1 #[1] shift(::Array{Any,1}, ::Array{Any,1}, ::Int64, ::Int64) at .\In[69]:11
            #motifs_B_i becomes postion + 1 until motif length
            if position+k-1 > length(dna[i])
                temp[i] = motifs_B[i]
            else
                temp[i] = SubString(dna[i],position,position+k-1)
                #change backgrounds
                if position == 1
                    prefix = ""
                    suffix = SubString(dna[i],position+k,length(dna[i]))
                elseif position == length(dna[i])-k+1
                    prefix = SubString(dna[i],1,position-1)
                    suffix = ""
                else
                    prefix = SubString(dna[i],1,position-1)
                    suffix = SubString(dna[i],position+k,length(dna[i]))
                end
                temp_backgrounds[i] = string(prefix,"",suffix)
            end
        end
        score2 = product(temp, backgrounds_B, t, k)
        if max(score1, score2) == score2
            motifs_B = temp
            backgrounds_B = temp_backgrounds
        end
    end
    for j in 1:5
        score1 = product(motifs_B, backgrounds_B, t, k)
        for i in 1:t
            #return 1:2
            temp = motifs_B
            temp_backgrounds = backgrounds_B
            position = findfirst(motifs_B[i], dna[i])[1] - 1
            #motifs_B_i becomes postion + 1 until motif length
            if position == 0
                temp[i] = motifs_B[i]
            else
                temp[i] = SubString(dna[i],position,position+k-1)
                #change backgrounds
                if position == 1
                    prefix = ""
                    suffix = SubString(dna[i],position-1+k,length(dna[i]))
                elseif position == length(dna[i])-k+1
                    prefix = SubString(dna[i],1,position-1)
                    suffix = ""
                else
                    prefix = SubString(dna[i],1,position-1)
                    suffix = SubString(dna[i],position+k,length(dna[i]))
                end
                temp_backgrounds[i] = string(prefix,"",suffix)
            end
        end
        score2 = product(temp, backgrounds_B, t, k)
        if max(score1, score2) == score2
            motifs_B = temp
            backgrounds_B = temp_backgrounds
        end
    end
    #println(motifs_B)
    #println(backgrounds_B)
    return(motifs_B, backgrounds_B)
end
#typeof(findfirst("ACT", "TACTAAAT")[1])

function Getmotifs(dna,k,t)
    motifs = []
    backgrounds = []
    for i in 1:t
        randomnumber = rand(1:length(dna[i])-k+1)
        tempmotif = SubString(dna[i],randomnumber,randomnumber+k-1)
        motifs = vcat(motifs, tempmotif)
        if randomnumber == 1

            prefix = ""
            suffix = SubString(dna[i],randomnumber+k,length(dna[i]))

        elseif randomnumber == length(dna[i])-k+1
            prefix = SubString(dna[i],1,randomnumber-1)
            suffix = ""

        else
            prefix = SubString(dna[i],1,randomnumber-1)
            suffix = SubString(dna[i],randomnumber+k,length(dna[i]))

        end
        background = string(prefix,"",suffix)
        backgrounds = vcat(backgrounds,background)
    end

    return(motifs,backgrounds)
end

function cross(t)
    crossoverpoint = rand(2:t-1)
    return(crossoverpoint)
end

function MotifFreq(motif, t, k)
    
    Acount = zeros(1,k)
    Ccount = zeros(1,k)
    Gcount = zeros(1,k)
    Tcount = zeros(1,k)
    for i in 1:k
        for j in 1:t
            if motifs[j][i] == 'A'
                Acount[i] = Acount[i] + 1
            elseif motifs[j][i] == 'C'
                Ccount[i] = Ccount[i] + 1
            elseif motifs[j][i] == 'G'
                Gcount[i] = Gcount[i] + 1
            else
                Tcount[i] = Tcount[i] + 1
            end
        end
        Acount[i] = Acount[i] + 1
        Ccount[i] = Ccount[i] + 1
        Gcount[i] = Gcount[i] + 1
        Tcount[i] = Tcount[i] + 1
    end
    
    Amotfreq = zeros(1,k)
    Cmotfreq = zeros(1,k)
    Gmotfreq = zeros(1,k)
    Tmotfreq = zeros(1,k)
    for i in 1:k
        totalcount = Acount[i] + Ccount[i] + Gcount[i] + Tcount[i]
        Amotfreq[i] = Acount[i]/totalcount
        Cmotfreq[i] = Ccount[i]/totalcount
        Gmotfreq[i] = Gcount[i]/totalcount
        Tmotfreq[i] = Tcount[i]/totalcount
    end
    return(Amotfreq,Cmotfreq,Gmotfreq,Tmotfreq)
    
end

function geneticAlgorithmMotifFinder(dna, k, t, G, r)
    population = 500
    generations = G
    motifMatrix = []
    backgroundMatrix = []
    for i in 1:population
        tempMotif = Getmotifs(dna,k,t)
        #println(tempMotif[1])
        #println(length(motifMatrix))
        if i == 1
            motifMatrix = tempMotif[1]
            #println(motifMatrix)
            backgroundMatrix = tempMotif[2]
        else
            motifMatrix = vcat(motifMatrix, tempMotif[1])
            backgroundMatrix = vcat(backgroundMatrix, tempMotif[2])
        end
    end

    #println("original " , motifMatrix)
    #check, seqPos, tempMotif, tempBack = mutation(motifMatrix[1:t])

    for generationCount in 1:generations
        #do the mutations at prob p for each set of motifs
        start = 1
        stop = t
        #mutations
        for i in 1:population
            #motif start - stop = new start - stop
            check, seqPos, tempMotif, tempBack = mutation(motifMatrix[start:stop], t, r, k, dna)
            if check == true
                motifMatrix[start+seqPos] = tempMotif
                backgroundMatrix[start+seqPos] = tempBack
            end
            start = start + t
            stop = stop + t
        end
        #MotifFreq(motifsA)[1]
        #motifc1 = motifsA[1:crossoverpoint] 
        #motifc1 = vcat(motifc1, motifsB[crossoverpoint+1:t])

        #pair up parents and make babies! children using cross overs
        parents = []
        babies = []
        babyBack = []
        for i in 1:(population/2)
            rand1 = rand(1:population)
            rand2 = rand(1:population)
            while rand1 in parents
                rand1 = rand(1:population)
            end
            while (rand1 == rand2) || (rand2 in parents)
                rand2 = rand(1:population)
            end
            push!(parents, rand1)
            push!(parents, rand2)
            start1 = ((rand1-1)*t)+1
            stop1 = start1 + t - 1
            start2 = ((rand2-1)*t)+1
            stop2 = start2 + t - 1

            crossoverpoint = cross(t)

            motifc1 = motifMatrix[start1:crossoverpoint+start1-1] 
            motifc1 = vcat(motifc1, motifMatrix[start2+crossoverpoint:stop2])

            motifc2 = motifMatrix[start2:crossoverpoint+start2-1] 
            motifc2 = vcat(motifc2, motifMatrix[start1+crossoverpoint:stop1])

            backc1 = backgroundMatrix[start1:crossoverpoint+start1-1] 
            backc1 = vcat(backc1, backgroundMatrix[start2+crossoverpoint:stop2])

            backc2 = backgroundMatrix[start2:crossoverpoint+start2-1] 
            backc2 = vcat(backc2, backgroundMatrix[start1+crossoverpoint:stop1])

            babies = vcat(babies, motifc1)
            babies = vcat(babies, motifc2)
            babyBack = vcat(babyBack, backc1)
            babyBack = vcat(babyBack, backc2)
        end
        #print("babies: ", babies)
        tournament = vcat(motifMatrix, babies)
        tournamentBack = vcat(backgroundMatrix, babyBack)
        #println(length(tournament))
        pairs = []
        productthingVector1 = []
        productthingVector2 = []
        newMatrix = []
        newBack = []
        #HUNGER GAMES
        for i in 1:population
            rand1 = rand(1:population*2)
            rand2 = rand(1:population*2)
            while rand1 in pairs
                rand1 = rand(1:population*2)
            end
            while (rand1 == rand2) || (rand2 in pairs)
                rand2 = rand(1:population*2)
            end
            start1 = ((rand1-1)*t)+1
            stop1 = start1 + t - 1
            start2 = ((rand2-1)*t)+1
            stop2 = start2 + t - 1
            productthingVector1 = product(tournament[start1:stop1], tournamentBack[start1:stop1], t, k)
            productthingVector2 = product(tournament[start2:stop2], tournamentBack[start2:stop2], t, k)
            #println("A score ", productthingVector1)
            #println("B score ", productthingVector2)
            #println("A ", tournament[start1:stop1])
            #    println("B ", tournament[start2:stop2])
            if productthingVector1 > productthingVector2
                newMatrix = vcat(newMatrix, tournament[start1:stop1])
                newBack = vcat(newBack, tournamentBack[start1:stop1])
                #println("A")
            else
                newMatrix = vcat(newMatrix, tournament[start2:stop2])
                newBack = vcat(newBack, tournamentBack[start2:stop2])
                #println("B")
            end
            #println("NEWMATRIX:    ", newMatrix)
        end
        motifMatrix = newMatrix
        backgroundMatrix = newBack
        #println(newMatrix)
        #get the probabilities 
        #vector = product(motifMatrix[1:t], backgroundMatrix[1:t])
        #println(vector)
        #println(productthingVector)
        #println(length(motifMatrix))
        #print(motifMatrix)
    end

    start1 = 1
    stop1 = start1 + t - 1
    start2 = t + 1
    stop2 = start2 + t - 1
    motifs_genConsensus = motifMatrix[start1:stop1]
    background_genConsensus = backgroundMatrix[start1:stop1]
    productthingMax = product(motifMatrix[start1:stop1], backgroundMatrix[start1:stop1], t, k)
    for i in 1:population-1
        productthingVector1 = product(motifMatrix[start1:stop1], backgroundMatrix[start1:stop1], t, k)
        productthingVector2 = product(motifMatrix[start2:stop2], backgroundMatrix[start2:stop2], t, k)

        if productthingVector1 > productthingMax
            motifs_genConsensus = motifMatrix[start1:stop1]
            background_genConsensus = backgroundMatrix[start1:stop1]
                #println("A")
        elseif productthingVector2 > productthingMax
            motifs_genConsensus = motifMatrix[start2:stop2]
            background_genConsensus = backgroundMatrix[start2:stop2]
        end
        start1 = start2
        stop1 = stop2
        start2 = stop1 + 1
        stop2 = start2 + t - 1
    end
    #println(motifs_genConsensus)

    global motifs_genAdj, background_genAdj = adjust(motifs_genConsensus, background_genConsensus, t, k, dna)
    #println(motifs_genAdj)

    motifs_Shift, background_Shift = shift(motifs_genAdj, background_genAdj, t, k, dna)
    #println(motifs_Shift)
    motifs_Consensus = motifs_Shift
    background_Consensus = background_Shift
    return(motifs_Consensus, background_Consensus)
end

function pwmScan(motifs_X, backgrounds_X, k, t, dna)
    motifs_A = deepcopy(motifs_X)
    backgrounds_A = deepcopy(backgrounds_X)
    global new_back = ""
    global output = ""
    max_prob = 0
    secondsites = []
    location = []
    for m in 1:t
        #println(m)
        DNA_Sequence = deepcopy(motifs_A[m])
        motifs_A_notI = deepcopy(motifs_A[1:end .!= m])
        backgrounds_A_notI = deepcopy(backgrounds_A[1:end .!= m])
        Acount, Ccount, Gcount, Tcount = bpCount(motifs_A_notI, t-1, k)
        Amotfreq, Cmotfreq, Gmotfreq, Tmotfreq = motFreq(Acount, Ccount, Gcount, Tcount, t-1, k)
        #println(Amotfreq)
        #println(Cmotfreq)
        #println(Gmotfreq)
        #println(Tmotfreq)
        prob_start = 1
        for j in 1:length(DNA_Sequence)
                if DNA_Sequence[j] == 'A'
                    prob_start = prob_start*Amotfreq[j]
                #println(prob_start)
                elseif DNA_Sequence[j] == 'C'
                    prob_start = prob_start*Cmotfreq[j]
                #println(prob_start)
                elseif DNA_Sequence[j] == 'G'
                    prob_start = prob_start*Gmotfreq[j]
                #println(prob_start)
                else
                    prob_start = prob_start*Tmotfreq[j]
                #println(prob_start)
                end
        end
        prob_start = prob_start * .9
        if m == 1
            max_prob = prob_start
                    #println(prob_start)
        end
        score1 = product(motifs_A_notI, backgrounds_A_notI, t-1, k)
        temp_backgrounds = ""
        subsetvec = []
        probvec = []
        global output = ""
        global new_back = ""
        #new_motif_location = A vector that is number of DNA sequences
        start_pos = findfirst(motifs_A[m], dna[m])[1]
        stop_pos = start_pos + k -1
        i = 1
        saved_pos = 0
        while i < length(dna[m])-k+1
            if ((i+k-1) < start_pos) || (i > stop_pos)
                
               # println("i ", i) #SubString(dna[i],i,i+k-1)
                #println("DNA[m]", dna[m])

                subsetDNA = SubString(dna[m],i,i+k-1)
                #println("subsetDNA  ", subsetDNA)
                prob = 1
                for j in 1:length(subsetDNA)
                    if subsetDNA[j] == 'A'
                        prob = prob*Amotfreq[j]
                    elseif subsetDNA[j] == 'C'
                        prob = prob*Cmotfreq[j]
                    elseif subsetDNA[j] == 'G'
                        prob = prob*Gmotfreq[j]
                    else
                        prob = prob*Tmotfreq[j]
                    end
                end
                    if i == 1
                        prefix = ""
                        suffix = SubString(dna[m],i+k,length(dna[m]))
                    elseif i == length(dna[m])-k+1
                        prefix = SubString(dna[m],1,i-1)
                        suffix = ""
                    else
                        prefix = SubString(dna[m],1,i-1)
                        suffix = SubString(dna[m],i+k,length(dna[m]))
                    end
                    temp_backgrounds = string(prefix,"",suffix)
                push!(subsetvec, subsetDNA)
                push!(probvec, prob)
                #println(prob)
                #println(maximum(probvec))
                if i == 1 && prob > max_prob
                    #println(subsetDNA)
                    global output = deepcopy(subsetDNA)
                    global new_back = temp_backgrounds
                    saved_pos = i
                    max_prob = prob
                    #println(output)
                else
                    if prob == maximum(probvec) && prob > max_prob
                        global output = deepcopy(subsetDNA)
                        global new_back = deepcopy(temp_backgrounds)
                        max_prob = prob
                        saved_pos = i
                        #println(prob)
                        #println(output)
                        #new_motif_location[random_number] = i                
                    end
                end
            end
            i = i + 1
        end
        #println(random)
        #println("-")
        
        if output != ""
            #println("output \t" ,output)
            #println("motif before whatever ",motifs_A[m])
            motifs_A[m] = deepcopy(output)
            #println("motif after whatever ",motifs_A[m])
            backgrounds_A[m] = new_back
            #println(motifs_A)
            #push!(secondsites, output)
            push!(location, saved_pos)
        else
            push!(location, start_pos)
        end


        ##un this a million times over

    end
    #println(motifs_A)
    new_motifs = deepcopy(motifs_A)
    #println(motifs_X)
    for i in 1:t
        if motifs_A[i] == motifs_X[i]
            if location[i] == findfirst(motifs_A[i], dna[i])[1]
                new_motifs[i] = "There is no second motif in this sequence"
            end
        end
    end
    
    return(motifs_A, new_motifs, location)
end

function callGame(dna, k, t, G, r, PWN)
    motifs_Consensus, background_Consensus = geneticAlgorithmMotifFinder(dna, k, t, G, r)
    #println("Monsensus motif\t",motifs_Consensus)
    #println(motif)
    println("These are the consensus motifs:")
    for i in 1:t
        println("For sequence ", i)
        println(motifs_Consensus[i], " at position ", findfirst(motifs_Consensus[i], dna[i])[1])
    end
    if PWN == true
        motif, new, loc = pwmScan(motifs_Consensus, background_Consensus, k, t, dna)
        println("")
        println("These are the secondary motifs:")
        for i in 1:t
            if new[i] != "There is no second motif in this sequence"
                println("For sequence ", i)
                println(new[i], " at position ", loc[i])
            end
        end
    end
    Acount, Ccount, Gcount, Tcount = bpCount(motifs_Consensus, t, k)
    Amotfreq, Cmotfreq, Gmotfreq, Tmotfreq = motFreq(Acount, Ccount, Gcount, Tcount, t, k)
    println("Probability vector for A: ", Amotfreq)
    println("Probability vector for C: ", Cmotfreq)
    println("Probability vector for G: ", Gmotfreq)
    println("Probability vector for T: ", Tmotfreq)

end

end # module
