import streamlit as st
import numpy as np
import pandas as pd

# Scoring constants for nucleotide sequences
MATCH = 1        # Match score for biological relevance
MISMATCH = -1    # Mismatch penalty
GAP_OPEN = -2    # Gap opening penalty
GAP_EXTEND = -2  # Gap extension penalty

# --- TABS ---
tabs = st.tabs(["🏠 HOME PAGE", "ℹ️ About", "👤 Team"])

# --- HOME PAGE ---
with tabs[0]:
    st.title("🧬 EvoAlign - Align, Visualise, Analyse")

    st.subheader("📥 Input FASTA Sequences")
    example_fasta = ">Seq1\nACTGGTAC\n>Seq2\nACTGTTAC\n>Seq3\nACTAGTAC\n>Seq4\nACTGGGAC"
    fasta_input = st.text_area("Enter DNA sequences in FASTA format", height=200, placeholder=example_fasta)

    def parse_fasta(text):
        sequences = {}
        label = None
        for line in text.strip().splitlines():
            if line.startswith(">"):
                label = line[1:]
                sequences[label] = ""
            else:
                sequences[label] += line.strip().upper()
        return sequences

    def get_score(a, b):
        if a == "-" or b == "-":
            return GAP_EXTEND
        return MATCH if a == b else MISMATCH

    def needleman_wunsch(seq1, seq2):
        m, n = len(seq1), len(seq2)
        score_mat = np.zeros((m + 1, n + 1), dtype=int)
        traceback = np.zeros((m + 1, n + 1), dtype=str)

        # Initialize with gap penalties
        for i in range(1, m + 1):
            score_mat[i][0] = GAP_OPEN + (i - 1) * GAP_EXTEND
            traceback[i][0] = "↑"
        for j in range(1, n + 1):
            score_mat[0][j] = GAP_OPEN + (j - 1) * GAP_EXTEND
            traceback[0][j] = "←"
        traceback[0][0] = "↖"

        # Fill the matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                diag = score_mat[i - 1][j - 1] + get_score(seq1[i - 1], seq2[j - 1])
                up = score_mat[i - 1][j] + (GAP_EXTEND if traceback[i - 1][j] == "↑" else GAP_OPEN)
                left = score_mat[i][j - 1] + (GAP_EXTEND if traceback[i][j - 1] == "←" else GAP_OPEN)
                score_mat[i][j] = max(diag, up, left)

                if score_mat[i][j] == diag:
                    traceback[i][j] = "↖"
                elif score_mat[i][j] == up:
                    traceback[i][j] = "↑"
                else:
                    traceback[i][j] = "←"

        # Traceback
        aligned1, aligned2 = "", ""
        i, j = m, n
        while i > 0 or j > 0:
            direction = traceback[i][j]
            if direction == "↖":
                aligned1 = seq1[i - 1] + aligned1
                aligned2 = seq2[j - 1] + aligned2
                i -= 1
                j -= 1
            elif direction == "↑":
                aligned1 = seq1[i - 1] + aligned1
                aligned2 = "-" + aligned2
                i -= 1
            elif direction == "←":
                aligned1 = "-" + aligned1
                aligned2 = seq2[j - 1] + aligned2
                j -= 1

        # Calculate final score from alignment
        final_score = 0
        gap_open = False
        for a, b in zip(aligned1, aligned2):
            if a == "-" or b == "-":
                if not gap_open:
                    final_score += GAP_OPEN
                    gap_open = True
                else:
                    final_score += GAP_EXTEND
            else:
                final_score += MATCH if a == b else MISMATCH
                gap_open = False

        return score_mat, final_score, aligned1, aligned2

    def format_alignment(a1, a2):
        match_line = "".join("|" if x == y else " " for x, y in zip(a1, a2))
        return f"Seq1: {a1}\n      {match_line}\nSeq2: {a2}"

    def calculate_distance_matrix(seqs, names):
        n = len(seqs)
        dist_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i, n):
                if i == j:
                    dist_matrix[i][j] = 0
                else:
                    _, _, a1, a2 = needleman_wunsch(seqs[i], seqs[j])

                    matches = sum(x == y for x, y in zip(a1, a2))
                    alignment_length = len(a1)
                    identity = (matches / alignment_length) * 100

                    dist_matrix[i][j] = dist_matrix[j][i] = 1 - (identity / 100)

        return dist_matrix

    if fasta_input:
        sequences = parse_fasta(fasta_input)
        names = list(sequences.keys())
        seqs = list(sequences.values())

        valid_chars = set("ACGT")
        invalid = any(c not in valid_chars for seq in seqs for c in seq)

        if invalid:
            st.error("❌ Invalid DNA characters detected. Allowed: A, C, G, T.")
        elif len(seqs) < 2:
            st.error("⚠️ Please enter at least 2 sequences.")
        else:
            # Calculate distance matrix
            dist_matrix = calculate_distance_matrix(seqs, names)

            st.subheader("📊 Pairwise Distance Matrix")
        
            dist_df = pd.DataFrame(dist_matrix, index=names, columns=names)
            st.dataframe(dist_df)

            st.subheader(f"🧬 Global Alignment Results")
            st.markdown(f"**Scoring Parameters:** Match = {MATCH}, Mismatch = {MISMATCH}, "
                       f"Gap Open = {GAP_OPEN}, Gap Extend = {GAP_EXTEND}")

            best_pair = ("", "")
            best_identity = -1
            worst_pair = ("", "")
            worst_identity = 101

            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    s1, s2 = seqs[i], seqs[j]
                    mat, score, a1, a2 = needleman_wunsch(s1, s2)

                    matches = sum(x == y for x, y in zip(a1, a2))
                    alignment_length = len(a1)
                    identity = (matches / alignment_length) * 100

                    if identity > best_identity:
                        best_identity = identity
                        best_pair = (names[i], names[j])
                    if identity < worst_identity:
                        worst_identity = identity
                        worst_pair = (names[i], names[j])

                    st.markdown(f"### 🔬 {names[i]} vs {names[j]}")
                    st.markdown(f"**Identity:** `{identity:.1f}%` ({matches}/{alignment_length} bases)")
                    st.code(format_alignment(a1, a2), language="text")

                    # Displaying alignment matrix
                    st.subheader(f"Alignment Matrix: {names[i]} vs {names[j]}")
                    st.write(pd.DataFrame(mat))

            st.subheader("🔍 Conclusion")
            col1, col2 = st.columns(2)
            with col1:
                st.success(f"**Most Similar Pair**\n\n"
                          f"**{best_pair[0]}** & **{best_pair[1]}**\n\n"
                          f"Identity: {best_identity:.1f}%")
            with col2:
                st.error(f"**Most Distant Pair**\n\n"
                         f"**{worst_pair[0]}** & **{worst_pair[1]}**\n\n"
                         f"Identity: {worst_identity:.1f}%")

# --- ABOUT PAGE ---
with tabs[1]:
    st.title("About This App")
    st.markdown("""    
    Tool Overview:
    EvoAlign is a bioinformatics tool designed to align DNA sequences and visualize their relationships through the global alignment algorithm. It provides an intuitive interface for analyzing nucleotide sequences in FASTA format, helping users understand genetic similarity through pairwise alignment and distance matrix generation.
    """)

    st.markdown("""
    **Objectives:**
    - To provide an easy-to-use interface for DNA sequence alignment and analysis.
    - To ensure accurate pairwise comparisons using the Needleman-Wunsch global alignment algorithm.
    - To facilitate the exploration of sequence relationships through dynamic distance matrix visualization.
    """)

    st.markdown("""
    **Features:**
    - Global alignment (Needleman-Wunsch)
    - Pairwise distance matrix
    - Detailed alignment visualization
    - Interactive matrix exploration
    """)

# --- TEAM PAGE ---
with tabs[2]:
    st.title("Meet the Developer")
    st.markdown("""
    Hi! I'm Shreeya,
    MSc Bioinformatics student passionate about genomics and neuroscience.

    **Acknowledgements:**
    I express my heartfelt gratitude to my mentors Dr. Kushagra Kashyap and Dr. Poonam Deshpande for their guidance and Streamlit for providing the platform to build this tool.
    """)
