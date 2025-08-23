---
title: 'An Introduction to BLAST (Basic Local Alignment Search Tool) for non-Bioinformaticians'
date: 2024-11-13
permalink: /posts/2024/11/blast/
tags:
  - bioinformatics
  - BLAST
  - science
---

Imagine you're a Geneticist at NASA, doing your rounds alone in the vast, red landscape, with only your sampling swabs for company. It's been hours already, and you're starting to get hungry, so you come back to your ship with your shoulders down, expecting more negative results and a disappointing dose of spaceship food. But when you test the last swab, you stumble across something extraordinary  – a sample that tests positive for DNA! Have you found life on Mars?

[Image by Nicolas Lobos (Unsplash)](images/post1-mars.png)

A wave of excitement rushes through you, and you run to your portable sequencer and after following the necessary protocols (forgetting your hunger!) you excitedly open the .FASTA file with the nucleotidic sequence of your sample. Beautiful patterns of A, C, G, and T illuminate your face through the screen:

    TCTTCGCAGGCTTCTTAATCTCAAACCTCATCCCGCCCCTACTAGTTCCACAAACAACCATACCAACTTACATAAAACTCGCAGCCCTCACCGTCACAGCAGCAGGATTCATTCTCGCCATAGAATTAAATCAAATCACA
    CTTCACATCAAAGACACCCGTCCCACACACATACACTCATTCTCATCCCTCTTAGGATACTACCCCAATGTCATGCACCGCCTGGCCCCCTTCCACACCCTTTCAATAAGCCAAAACCTGGCATCCCTCTTAGACCTACT
    ATGACTAGAAAAAGCCATTCCCAAAAGCCTTTCTCAACTTCAACTTTCAGCCTCCGCAACCACATCCAACCAAAAAGGACTAATTAAACTCTACTTCTTATCTTTTCTCCTCTCCCTCTCACTAGGCCTACTAATCCTCC
    TCTAACGCCCACGAGTAATTTCGATCACAATAAAAATACTCACAAATAAAGATCACCCGGCCACAACCACCAATCAGCTACCATAACTATACAAGGCGGATCCTCCAATATAATCCTCACGAACCAAACTCATATCATCA
    CCCCCTAAAATAACCCAGTCCCCTATATCCTTTAAACCACAAAACACAAGCCCTATACCACTCTCACTTATCAGTCAAACTACCATCAACGCCTCCGCCAATAAGCCTACAATCAAACCCCCCAAAATCACACTATCTGA
    CCCCCAAGCTTCAGGATATTCTTCAGTAGCCATGGCGGTAGTGTAACCAAACACTACCAACATGCCCCCCAGATAAACCAAAAATACTATCAACCCCAGAAAAGACCCGCCGAGACTCATTACAATACCACACCCAACTC
    CACCGCCCACAATTAACCCCAACCCACCATATACAGGGGATGGCTTTGAAGAAAAACTAACAAAGCTGACCACAAGGATAACACTCAATAAAAACACCATATACGTCATAGTTCCCGCATGGACCTAACCATGACCAATG
    ATATGAAAAACCACCGTTGTAATTTCAACTACAAGAACCAATGACCAACATCCGCAAAACTCATCCTCTCTTTAAAATTATCAACCACTCATTTATTGACCTCCCCACCCCAACAAGTATTTCAGCATGATGGAACTTCG
    GCTCCCTACTAGGCATTTGCCTTCTCATCCAAATCGTCACAGGCCTATTTCTAGCAATACACTATACATCAGACACACTCACCGCCTTCTCATCCGTTACCCACATTTGCCGAGACGTAAATTACGGATGAATCATCCGC
    TACATGCATGCCAACGGAGCGTCTCTATTCTTCATATGCCTTTACCTCCACGTAGGCCGTGGAATGTACTATGGGTCCTACACATTCACGGAAACATGAAACATCGGAGTAGTACTTCTACTAACAGTCATAGCCACAGC
    ATTCATGGGATACGTCCTCCCATGAGGGCAAATATCTTTCTGAGGAGCTACAGTAATCACCAACCTCCTATCCGCTATTCCATACATCGGAACTGACCTAGTCGAGTGAATCTGAGGTGGGTTTTCAGTAGACAAAGCGA
    CCTTAACACGATTCTTCGCCTTCCACTTTATCCTCCCCTTTGTCGTTACAGCCCTAGTGATAATTCACCTACTGTTCCTACACGAAACAGGGTCCAACAACCCAACTGGCATATCATCTACTATAGATGCAATCCCATTT
    CACCCGTACTACACTATTAAAGACATTCTAGGCCTATTTCTCATGATCCTCTTCCTAATAACGCTAGTTCTATTTGCCCCAGATCTCCTAGGGGACCCAGATAACTACATCCCAGCAAACCCATTAAGCACACCTCCCCA

But like this, the sequence is almost meaningless to you. Indecipherable. What if your sampled Martian is related to some species back home on Earth? Or... what if it's actually just contamination from something you brought with you? You need to find out.

Being a Geneticist, you know that "homology" in biology refers to the degree of similarity between two things due to a common ancestry. For example, the human tailbone is homologous to the tails of other primates, because they have a common origin. As with anatomical structures, sequence similarity between two proteic or DNA sequences may imply a common origin, which can occur in different species (in which case we say that both sequences are orthologs of each other), or in the same species by a duplication event (in which case the repeated sequences are called paralogs). This, in principle, is the main idea behind DNA paternity tests. Your DNA is similar to that of your biological father.

All these things considered, you go to your database of all the genomes sequenced on Planet Earth, and open each file one by one to check for matches with your sampled sequence...

...not. You run [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), the algorithm that could help you solve this puzzle. It is a widely used tool for finding regions of similarity between biological sequences. By comparing your mysterious sample against known genetic sequences, it can help you uncover whether this Martian organism shares a hidden connection with life on our home planet.

Alright, that's enough introduction. But how does the algorithm work? Let's go over it in **four steps: indexing, local alignment, extension, and statistical significance filtering.**

Step 1:  Indexing
======

Our query sequence (Q) is broken down, by sliding window, into smaller sub-sequences that we call "words" or k-mers, of length k. Each word is then assigned an index according to its position in Q, and a unique numerical key. For each word, BLAST uses a numerical scoring system called a [scoring/substitution matrix](https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/2105/index.php?manual=BE_Scoring_matrices.html) (such as [BLOSUM62](https://en.wikipedia.org/wiki/BLOSUM)) to generate a list of "neighboring words" that are similar to it.

Example: For a DNA sequence like "AGACTTCA," let's examine a single-base substitution, focusing on the first 'A' (adenine). Adenine, a two-ring structure classified as a purine, can be replaced by either the other purine, guanine (G), or by a one-ring pyrimidine, such as cytosine (C) or thymine (T). Substitutions between bases with similar structures, such as purine-to-purine or pyrimidine-to-pyrimidine, are called transitions. Conversely, substitutions between bases with differing ring structures, such as purine-to-pyrimidine or vice versa, are known as transversions. Substitution matrices penalize transversions more heavily than transitions because the structural and chemical differences are more pronounced, potentially leading to greater functional changes in the sequence.

[Transitions and transversions: the DNA's mutations](images/post1-transitions_transversions.png)

The simplest substitution matrix I can think of is one where a transition deducts 1 point, a transversion deducts 2 points, and no change in the base adds 1 point. For example, if the query sequence " AGACTTCA" aligns with the reference sequence "GGACTTCA", the score would be calculated as: -1 (for the transition) + 1 + 1 + 1 + 1 + 1 + 1 + 1 = 6 points. If the sequence changes to "TGACTTCA", the score would be: -2 (for the transversion) + 1 + 1 + 1 + 1 + 1 + 1 + 1 = 5 points. This demonstrates the basic concept behind substitution scoring matrices.

[My simple DNA substitution matrix where the column nucleotides originate from the query sequence, and the row nucleotides are the possible substitutions.](images/post1-matrix.png)

A neighboring word is only included in the list if the similarity score (based on the chosen scoring matrix) exceeds a certain predetermined **threshold_1** (remember this parameter).These high-scoring words are stored in a data structure called a [hash-table](https://en.wikipedia.org/wiki/Hash_table), which allows for fast lookups similarly to how a telephone directory would. The indexing is precomputed for sequences that are already uploaded in a database, to speed up the process when running BLAST.

By following these steps, BLAST efficiently identifies and ranks local alignments between a query sequence and database sequences. Its design balances speed and sensitivity, making it highly effective for large-scale sequence searches.

[example table with k-words, their position, and their numerical key](images/post1-kmers.png)

Step 2: Local Alignment
======

This step is all about finding hits, and it's where the "LA" part of BLAST comes in. A local alignment, in contrast to a global one, is an algorithm that tries to find the best matching sub-sequences between two input sequences (here: our input Martian sample and each sequence in the general database of Earth's species), even if the full sequences aren't completely aligned. A similar type of local alignment method is the dynamic programming [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) algorithm. Being that such dynamic algorithms are extremely demanding computationally, it is very convenient that we have compressed our data in the previous step!

BLAST checks each position in your query sequence (Q[i], where i is an index that represents a particular position in the query sequence Q) for matches with the database. Exact matches between the database and lists of Q[i] are searched for and subsequently saved as seed hits. Here, I represented the different Q[i]s our sample sequence was broken down into, matching exactly with different regions of the full genomes of the database.

[A local alignment visualization in BLAST](images/post1-align.png)

Step 3: Extension
======

For each seed hit found, BLAST attempts to extend the alignment in both directions (left and right) to find a local alignment. The difference with the previous step, is that these alignments aren't necessarily exact, but allow for differences between the sequences. The algorithm continues aligning bases (for DNA) or amino acids (for proteins) until the score drops below a cutoff value that we can call **threshold_2** (remember this as well). The result of this step are the high-scoring sequence pairs (HSPs), which are gapped or ungapped alignments of significant length, i.e. the highest scoring points, where either shortening or lengthening the alignment would result in a worse sum total.

[Algorithm details](images/post1-hsp.png)

Step 4: Statistical significance filtering
======

After finding the HSPs, BLAST calculates their statistical significance to determine whether they are likely to have occurred by chance. This is done by computing the E-value (expected value), which indicates the number of alignments with an equal or better score expected by random chance. This computation takes into account: (a) the size of the database used, (b) the scoring system, and (c) the length of the query. This normalizes the score to make it possible to compare results, even when different substitution matrices were used.

Alignments with a low E-value (e.g., < 0.001) are considered statistically significant and reported as hits. All significant alignments are ranked and BLAST outputs the results, including information like alignment score, query coverage (the percentage the query sequence that is included in the reference sequence), and sequence identity (the percentage of exact coincidences between both sequences).

That's nice and all, but how come BLAST is so fast if so many steps are involved? How can it search through so many genomes at the same time and return an answer in a few seconds? The answer is simple: heuristics. Heuristic methods are proceedings that take pragmatic or informal steps to obtain a result faster, even if such answer might end up being sub-optimal. When we divided our sequence into "k-mers" ("words"), we simplified our starting sequence in order to reach an approximation of the best alignment. Nothing guarantees, like in the case of the Smith-Waterman algorithm I told you about before, that we reach the optimal answer. This makes BLAST a very fast and specific (the ability to correctly identify a true matching sequence is high), but the downside is that we lose specificity (we might find some false positives as well), and this is a limitation that you have to take into account. In practice, though, this algorithm works pretty neatly. There's a reason why upon publication, its article is [one of the most cited Papers to date](https://pubmed.ncbi.nlm.nih.gov/2231712/).

Another heuristic approximation that BLAST makes is assuming that if two sequences are similar, they have a "word" in common. This assumption helps the program run faster, but may not reflect a biological reality. A larger word size would mean a faster run, but it also becomes less sensitive to real results. When it comes to the threshold_1 and threshold_2 I mentioned before, the higher they are, the more stringent this method becomes, and so if we choose high values for them, BLAST would be slower, yet more sensitive. This makes sense, because if you are stricter, the chances of getting a false positive really reduces. In turn, the chances of false negatives greatly increase - it's all a game of balance.

I hope you found this explanation helpful! I’d love to hear your thoughts, questions, or suggestions—feel free to share them.

P.S.: Are you up for a challenge? Then use BLAST to analyze the sampled DNA sequence. Head to [NCBI's online tool](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and perform a "nucleotide-to-nucleotide" search, comparing your query sequence to matched DNA sequences. Explore which species the sequence is most closely related to, evaluate whether the statistics support that conclusion or not, and share your insights if you want to.
