#############################################################################
# Makefile for building: pagan2
#############################################################################

####### Compiler, tools and options
UNAME    = $(shell uname)
PWD      = /pagan

CC       = /usr/bin/gcc-4.8
CXX      = /usr/bin/g++-4.8
LEX      = flex
YACC     = yacc
CFLAGS   = -pipe -Wall -Wno-unused-parameter -Wno-unused -O2 -fopenmp 
CXXFLAGS = -std=c++11 -pipe -Wall -Wno-unused-parameter -Wno-unused -Wno-deprecated -O2 -fopenmp -DNCBI_TOOLKIT 
LEXFLAGS =
YACCFLAGS= -d
NCBIPATH = ../ncbi/GCC480-ReleaseMT64
NCBIINC  = -I$(NCBIPATH)/inc -I$(NCBIPATH)/../include
INCPATH  = -I../boost/include -I../curl/include -I. -I./utils -I./main -I/usr/include $(NCBIINC)
LINK     = /usr/bin/g++-4.8
LFLAGS   = -g -static 
NCBILIBS = -L../boost/lib -L../curl/lib -L$(NCBIPATH)/lib -lxblast-static -lxalgoblastdbindex-static -lcomposition_adjustment-static -lxalgodustmask-static -lxalgowinmask-static -lseqdb-static -lxobjutil-static -lxnetblast-static -lblastdb-static -lscoremat-static -ltables-static -lxconnect-static -lxobjmgr-static -lgenome_collection-static -lseqset-static -lseq-static -lseqcode-static -lsequtil-static -lpub-static -lmedline-static -lbiblio-static -lgeneral-static -lxser-static -lxutil-static -lxncbi-static -lpthread -ldl
LIBS     = $(SUBLIBS)  -Wl,-rpath, -lboost_program_options -lboost_regex -lboost_thread -lboost_system -lgomp -lm ../curl/lib/libcurl.a -lz $(NCBILIBS)
ifneq ($(UNAME), Darwin)
LIBS    += -lrt
endif
AR            = ar cqs
RANLIB        = 
TAR           = tar -cf
GZIP     = gzip -9f
COPY          = cp -f
COPY_FILE= $(COPY) -p
COPY_DIR = $(COPY) -pR
DEL_FILE      = rm -f
SYMLINK  = ln -sf
DEL_DIR       = rmdir
MOVE          = mv -f

####### Output directory

OBJECTS_DIR   = ./

####### Files

HEADERS       = main/node.h \
		main/sequence.h \
		main/basic_alignment.h \
		main/viterbi_alignment.h \
		main/reference_alignment.h \
		main/reads_aligner.h\
		utils/text_utils.h \
		utils/settings.h \
		utils/newick_reader.h \
		utils/model_factory.h \
		utils/fasta_reader.h \
		utils/eigen.h \
		utils/exceptions.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		utils/fasta_entry.h \
		utils/settings_handle.h \
		utils/xml_writer.h \
		utils/evol_model.h \
		utils/check_version.h \
		utils/log_output.h \
		utils/find_anchors.h \
		utils/codon_translation.h \
		utils/input_output_parser.h \
		utils/mafft_alignment.h \
		utils/raxml_tree.h \
		utils/tree_node.h \
		utils/bppdist_tree.h \
		utils/bppphysamp_tree.h \
		utils/bppancestors.h \
		utils/fasttree_tree.h \
		utils/ncbi_blast.h \
		utils/tunnel_matrix.h
SOURCES       = main.cpp \
		main/node.cpp \
		main/sequence.cpp \
		main/reference_alignment.cpp \
		main/viterbi_alignment.cpp \
		main/basic_alignment.cpp \
		main/reads_aligner.cpp \
		utils/text_utils.cpp \
		utils/settings.cpp \
		utils/newick_reader.cpp \
		utils/model_factory.cpp \
		utils/fasta_reader.cpp \
		utils/eigen.cpp \
		utils/db_matrix.cpp \
		utils/int_matrix.cpp \
		utils/settings_handle.cpp \
		utils/xml_writer.cpp \
		utils/evol_model.cpp \
		utils/check_version.cpp \
		utils/exonerate_queries.cpp \
		utils/log_output.cpp \
		utils/find_anchors.cpp \
		utils/codon_translation.cpp \
		utils/input_output_parser.cpp \
		utils/mafft_alignment.cpp \
		utils/raxml_tree.cpp \
		utils/tree_node.cpp \
		utils/bppdist_tree.cpp \
		utils/bppphysamp_tree.cpp \
		utils/fasttree_tree.cpp \
		utils/bppancestors.cpp \
		utils/ncbi_blast.cpp \
		utils/tunnel_matrix.cpp
OBJECTS       = main.o \
		node.o \
		sequence.o \
		reference_alignment.o \
		viterbi_alignment.o \
		basic_alignment.o \
		reads_aligner.o \
		text_utils.o \
		settings.o \
		newick_reader.o \
		model_factory.o \
		fasta_reader.o \
		eigen.o \
		db_matrix.o \
		int_matrix.o \
		settings_handle.o \
		xml_writer.o \
		evol_model.o \
		check_version.o \
		exonerate_queries.o \
		log_output.o \
		find_anchors.o \
		codon_translation.o \
		input_output_parser.o \
		mafft_alignment.o \
		raxml_tree.o \
		tree_node.o \
		bppdist_tree.o \
		bppphysamp_tree.o \
		bppancestors.o \
		fasttree_tree.o \
		ncbi_blast.o \
		tunnel_matrix.o
TARGET        = pagan2

first: all
####### Implicit rules

.SUFFIXES: .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) $(TARGET)


FORCE:

####### Compile

main.o: main.cpp utils/settings.h \
		utils/settings_handle.h \
		utils/log_output.h \
		utils/input_output_parser.h \
		utils/fasta_entry.h \
		utils/fasta_reader.h \
		utils/exceptions.h \
		main/node.h \
		main/sequence.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h \
		utils/codon_translation.h \
		main/reads_aligner.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

node.o: main/node.cpp main/node.h \
		utils/exceptions.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o node.o main/node.cpp

sequence.o: main/sequence.cpp utils/settings.h \
		utils/model_factory.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/evol_model.h \
		utils/settings_handle.h \
		utils/text_utils.h \
		utils/exceptions.h \
		main/sequence.h \
		utils/fasta_entry.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o sequence.o main/sequence.cpp

reference_alignment.o: main/reference_alignment.cpp main/reference_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		utils/exceptions.h \
		main/node.h \
		main/viterbi_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o reference_alignment.o main/reference_alignment.cpp

viterbi_alignment.o: main/viterbi_alignment.cpp main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		utils/exceptions.h \
		utils/codon_translation.h \
		utils/find_anchors.h \
		utils/substring_hit.h \
		utils/exonerate_queries.h \
		main/node.h \
		main/reference_alignment.h \
		utils/ncbi_blast.h \
		utils/tunnel_matrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o viterbi_alignment.o main/viterbi_alignment.cpp

basic_alignment.o: main/basic_alignment.cpp main/basic_alignment.h \
		utils/model_factory.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/evol_model.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		utils/exceptions.h \
		main/node.h \
		main/viterbi_alignment.h \
		main/reference_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o basic_alignment.o main/basic_alignment.cpp

reads_aligner.o: main/reads_aligner.cpp main/reads_aligner.h \
		utils/settings.h \
		utils/settings_handle.h \
		utils/model_factory.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/evol_model.h \
		utils/fasta_entry.h \
		utils/fasta_reader.h \
		utils/exceptions.h \
		main/node.h \
		main/sequence.h \
		main/viterbi_alignment.h \
		main/basic_alignment.h \
		main/reference_alignment.h \
		utils/exonerate_queries.h \
		utils/substring_hit.h \
		utils/text_utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o reads_aligner.o main/reads_aligner.cpp

text_utils.o: utils/text_utils.cpp utils/text_utils.h \
		utils/exceptions.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o text_utils.o utils/text_utils.cpp

settings.o: utils/settings.cpp utils/settings.h \
		utils/check_version.h \
		utils/log_output.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o settings.o utils/settings.cpp

newick_reader.o: utils/newick_reader.cpp utils/newick_reader.h \
		utils/exceptions.h \
		main/node.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h \
		utils/text_utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o newick_reader.o utils/newick_reader.cpp

model_factory.o: utils/model_factory.cpp utils/model_factory.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/evol_model.h \
		utils/settings.h \
		utils/settings_handle.h \
		utils/eigen.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o model_factory.o utils/model_factory.cpp

fasta_reader.o: utils/fasta_reader.cpp utils/fasta_reader.h \
		utils/exceptions.h \
		main/node.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h \
		utils/text_utils.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fasta_reader.o utils/fasta_reader.cpp

eigen.o: utils/eigen.cpp utils/eigen.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o eigen.o utils/eigen.cpp

db_matrix.o: utils/db_matrix.cpp utils/db_matrix.h \
		utils/log_output.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o db_matrix.o utils/db_matrix.cpp

int_matrix.o: utils/int_matrix.cpp utils/int_matrix.h \
		utils/log_output.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o int_matrix.o utils/int_matrix.cpp

settings_handle.o: utils/settings_handle.cpp utils/settings_handle.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o settings_handle.o utils/settings_handle.cpp

xml_writer.o: utils/xml_writer.cpp utils/xml_writer.h \
		utils/exceptions.h \
		main/node.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		utils/fasta_entry.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o xml_writer.o utils/xml_writer.cpp

evol_model.o: utils/evol_model.cpp utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/model_factory.h \
		utils/settings.h \
		utils/settings_handle.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o evol_model.o utils/evol_model.cpp

check_version.o: utils/check_version.cpp utils/check_version.h \
		utils/log_output.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o check_version.o utils/check_version.cpp

exonerate_queries.o: utils/exonerate_queries.cpp utils/exonerate_queries.h \
		utils/fasta_entry.h \
		utils/substring_hit.h \
		main/node.h \
		utils/exceptions.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o exonerate_queries.o utils/exonerate_queries.cpp

log_output.o: utils/log_output.cpp utils/log_output.h \
		utils/settings_handle.h \
		utils/settings.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o log_output.o utils/log_output.cpp

find_anchors.o: utils/find_anchors.cpp utils/find_anchors.h \
		utils/substring_hit.h \
		utils/settings_handle.h \
		utils/settings.h \
		main/sequence.h \
		utils/log_output.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o find_anchors.o utils/find_anchors.cpp

codon_translation.o: utils/codon_translation.cpp utils/codon_translation.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o codon_translation.o utils/codon_translation.cpp

input_output_parser.o: utils/input_output_parser.cpp utils/input_output_parser.h \
		utils/fasta_entry.h \
		utils/fasta_reader.h \
		utils/exceptions.h \
		main/node.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h \
		utils/codon_translation.h \
		utils/newick_reader.h \
		utils/mafft_alignment.h \
		utils/bppancestors.h \
		utils/bppdist_tree.h \
		utils/bppphysamp_tree.h \
		utils/raxml_tree.h \
		utils/exonerate_queries.h \
		utils/substring_hit.h \
		utils/xml_writer.h \
		utils/tree_node.h \
		main/reads_aligner.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o input_output_parser.o utils/input_output_parser.cpp

mafft_alignment.o: utils/mafft_alignment.cpp utils/mafft_alignment.h \
		utils/settings_handle.h \
		utils/settings.h \
		utils/fasta_entry.h \
		utils/log_output.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o mafft_alignment.o utils/mafft_alignment.cpp

raxml_tree.o: utils/raxml_tree.cpp utils/raxml_tree.h \
		utils/settings_handle.h \
		utils/settings.h \
		utils/fasta_entry.h \
		utils/model_factory.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/evol_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o raxml_tree.o utils/raxml_tree.cpp

tree_node.o: utils/tree_node.cpp utils/tree_node.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o tree_node.o utils/tree_node.cpp

bppdist_tree.o: utils/bppdist_tree.cpp utils/bppdist_tree.h \
		utils/settings_handle.h \
		utils/settings.h \
		utils/fasta_entry.h \
		utils/model_factory.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/evol_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o bppdist_tree.o utils/bppdist_tree.cpp

bppphysamp_tree.o: utils/bppphysamp_tree.cpp utils/bppphysamp_tree.h \
		utils/settings_handle.h \
		utils/settings.h \
		utils/fasta_entry.h \
		utils/log_output.h \
		utils/fasta_reader.h \
		utils/exceptions.h \
		main/node.h \
		main/sequence.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h \
		utils/newick_reader.h \
		utils/tree_node.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o bppphysamp_tree.o utils/bppphysamp_tree.cpp

bppancestors.o: utils/bppancestors.cpp utils/bppancestors.h \
		utils/fasta_entry.h \
		main/node.h \
		utils/exceptions.h \
		utils/settings.h \
		utils/settings_handle.h \
		main/sequence.h \
		main/viterbi_alignment.h \
		utils/evol_model.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		main/basic_alignment.h \
		utils/model_factory.h \
		main/reference_alignment.h \
		utils/fasta_reader.h \
		utils/codon_translation.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o bppancestors.o utils/bppancestors.cpp

fasttree_tree.o: utils/fasttree_tree.cpp utils/fasttree_tree.h \
		utils/settings_handle.h \
		utils/settings.h \
		utils/fasta_entry.h \
		utils/model_factory.h \
		utils/db_matrix.h \
		utils/log_output.h \
		utils/int_matrix.h \
		utils/evol_model.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fasttree_tree.o utils/fasttree_tree.cpp

ncbi_blast.o: utils/ncbi_blast.cpp utils/ncbi_blast.h \
		utils/log_output.h \
		utils/substring_hit.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ncbi_blast.o utils/ncbi_blast.cpp

tunnel_matrix.o: utils/tunnel_matrix.cpp utils/tunnel_matrix.h \
		main/basic_alignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o tunnel_matrix.o utils/tunnel_matrix.cpp

####### Install

install: all 

uninstall: 

