check:
	../examples/newick-fasta-unrooted $(abs_top_srcdir)/data/small.tre  $(abs_top_srcdir)/data/small.fas > nfu-output.txt && python $(abs_srcdir)/demand-all-lines.py $(abs_srcdir)/nfu-template.txt nfu-output.txt
	../examples/newick-fasta-rooted $(abs_top_srcdir)/data/small-rooted.tre  $(abs_top_srcdir)/data/small.fas > nfr-output.txt && python $(abs_srcdir)/demand-all-lines.py $(abs_srcdir)/nfr-template.txt nfr-output.txt
	../examples/partial $(abs_top_srcdir)/data/small.tre  $(abs_top_srcdir)/data/small.fas > partial-output.txt && python $(abs_srcdir)/demand-all-lines.py $(abs_srcdir)/partial-template.txt partial-output.txt
	../examples/protein-list $(abs_top_srcdir)/data/20.tre  $(abs_top_srcdir)/data/20.fas > protein-list.txt && python $(abs_srcdir)/demand-all-lines.py $(abs_srcdir)/protein-list-template.txt protein-list.txt

