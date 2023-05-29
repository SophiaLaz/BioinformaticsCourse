import random

def randomized_motif_search(Dna, k, t):
    def profile_most_probable_kmer(text, k, profile):
        n = len(text)
        max_prob = -1
        most_probable = ""
        for i in range(n - k + 1):
            kmer = text[i:i+k]
            prob = 1
            for j, nucleotide in enumerate(kmer):
                prob *= profile[nucleotide][j]
            if prob > max_prob:
                max_prob = prob
                most_probable = kmer
        return most_probable

    def create_profile(motifs, pseudocounts=False):
        profile = {'A': [], 'C': [], 'G': [], 'T': []}
        t = len(motifs)
        k = len(motifs[0])
        for nucleotide in ['A', 'C', 'G', 'T']:
            for j in range(k):
                count = sum((1 + pseudocounts) for motif in motifs if motif[j] == nucleotide)
                profile[nucleotide].append(count / (t + pseudocounts * 4))
        return profile

    def score(motifs):
        t = len(motifs)
        k = len(motifs[0])
        score = 0
        for j in range(k):
            count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            for i in range(t):
                count[motifs[i][j]] += 1
            max_count = max(count.values())
            score += t - max_count
        return score

    def randomized_motif_search_single(Dna, k, t):
        motifs = []
        for i in range(t):
            start = random.randint(0, len(Dna[i]) - k)
            motifs.append(Dna[i][start:start+k])

        best_motifs = motifs
        while True:
            profile = create_profile(motifs, pseudocounts=True)
            motifs = [profile_most_probable_kmer(seq, k, profile) for seq in Dna]
            if score(motifs) < score(best_motifs):
                best_motifs = motifs
            else:
                return best_motifs

    best_motifs = None
    best_score = float('inf')
    for _ in range(1000):
        motifs = randomized_motif_search_single(Dna, k, t)
        current_score = score(motifs)
        if current_score < best_score:
            best_motifs = motifs
            best_score = current_score

    return best_motifs

# Example usage
k = 8
t = 5
Dna = [
    "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]

best_motifs = randomized_motif_search(Dna, k, t)
for motif in best_motifs:
    print(motif)
