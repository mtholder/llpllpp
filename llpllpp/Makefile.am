pkglib_LTLIBRARIES = libllpllpp.la
libllpllpp_ladir = $(includedir)/llpllpp
libllpllpp_la_LDFLAGS = -release @PACKAGE_VERSION@ -lm -lpll

libllpllpp_la_HEADERS = \
	instance.hpp \
	tree.hpp

libllpllpp_la_SOURCES = \
	dsct_prob_model.cpp \
	instance.cpp \
	parsed_matrix.cpp \
	partitioned_data.cpp \
	phylo_calculator.cpp \
	tree.cpp