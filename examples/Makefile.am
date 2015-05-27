LDADD       = @top_builddir@/llpllpp/libllpllpp.la
AM_CPPFLAGS = -I@top_srcdir@/llpllpp
noinst_PROGRAMS = newick-fasta-unrooted \
	protein-list

newick_fasta_unrooted_SOURCES = newick-fasta-unrooted.cpp
newick_fasta_unrooted_CPPFLAGS = $(AM_CPPFLAGS)

protein_list_SOURCES = protein-list.cpp
protein_list_unrooted_CPPFLAGS = $(AM_CPPFLAGS)

check:
	echo 'no checks implemented yet'
