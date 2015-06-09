LDADD       = @top_builddir@/llpllpp/libllpllpp.la
AM_CPPFLAGS = -I@top_srcdir@/llpllpp
noinst_PROGRAMS = newick-fasta-optimize-one-edge \
	newick-fasta-unrooted \
	newick-fasta-rooted \
	partial \
	protein-list

partial_SOURCES = partial-traversal.cpp
partial_CPPFLAGS = $(AM_CPPFLAGS)

newick_fasta_optimize_one_edge_SOURCES = newick-fasta-optimize-one-edge.cpp
newick_fasta_optimize_one_edge_CPPFLAGS = $(AM_CPPFLAGS)

newick_fasta_rooted_SOURCES = newick-fasta-rooted.cpp
newick_fasta_rooted_CPPFLAGS = $(AM_CPPFLAGS)

newick_fasta_unrooted_SOURCES = newick-fasta-unrooted.cpp
newick_fasta_unrooted_CPPFLAGS = $(AM_CPPFLAGS)

protein_list_SOURCES = protein-list.cpp
protein_list_unrooted_CPPFLAGS = $(AM_CPPFLAGS)

check:
	echo 'no checks implemented yet'
