/***************************************************************************
 *   Copyright (C) 2010-2019 by Ari Loytynoja                              *
 *   ari.loytynoja@gmail.com                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include "utils/settings.h"
#include "utils/check_version.h"
#include "utils/log_output.h"

namespace po = boost::program_options;

using namespace std;
using namespace ppa;

Settings::Settings(){}

int Settings::read_command_line_arguments(int argc, char *argv[])
{
    version = 1.53;
    date = "29 August, 2019";

    boost::program_options::options_description minimal("Minimal progressive alignment options",100);
    minimal.add_options()
        ("seqfile,s", po::value<string>(), "sequence infile (FASTA)")
        ("treefile,t", po::value<string>(), "tree file")
    ;

    boost::program_options::options_description help_update("Help and updates",100);
    help_update.add_options()
        ("help,h", "display all program options")
        ("version,v","show program version and check for updates")
    ;

    boost::program_options::options_description generic("Generic options",100);
    generic.add_options()
        ("outfile,o", po::value<string>(), "alignment outfile")
        ("outformat,f", po::value<string>(), "alignment format [fasta,nexus,paml,phylipi,phylips,raxml]")
        ("xml,x","output XML alignment")
        ("events","output inferred evolutionary events")
        ("guidetree", "output alignment guidetree (with NHX tags)")
        ("ancestors", "include ancestors in outfile")
        ("translate", "translate DNA input to protein")
        ("mt-translate", "translate mtDNA input to protein")
        ("no-terminal-edges", "assume terminal gaps as missing data")
        ("silent","minimal output")
        ("config-file",po::value<string>(),"config file with additional arguments")
        ("config-log-file",po::value<string>(),"log file for given arguments")
        ("threads",po::value<int>(),"number of threads")
#ifdef NCBI_TOOLKIT
        ("no-ncbi", "do not use NCBI_toolkit")
#endif
    ;

    boost::program_options::options_description generic2("More generic options",100);
    generic2.add_options()
        ("xml-nhx","output XML alignment with NHX tree")
        ("raxml-tree","use RAxML for guidetree computation [default FastTree]")
        ("bppdist-tree","use BppDist for guidetree computation [default FastTree]")
        ("no-bppancestors","no BppAncestors (slow for large alignments)")
        ("noise", po::value<int>(), "output noise level")
        ("log-output-file",po::value<string>(),"output to file instead of stdout")
        ("temp-folder",po::value<string>(),"non-standard place for temp files")
    ;

    boost::program_options::options_description reads_alignment("Alignment extension options",100);
    reads_alignment.add_options()
        ("ref-seqfile,a", po::value<string>(), "reference alignment file (FASTA)")
        ("ref-treefile,r", po::value<string>(), "reference tree file (NH/NHX)")
        ("queryfile,q", po::value<string>(), "query file (FASTA/FASTQ)")
        ("fragments", "short queries: place together")
        ("both-strands","consider both strands, keep better (DNA)")
        ("find-orfs", "find ORFs, keep good (DNA)")
        ("fast-placement","use Exonerate to quickly assign queries to nodes")
        ("very-fast-placement","shorthand for very fast heuristic settings")
        ("all-nodes","place query to any node (default)")
        ("internal-nodes","place query to an internal node")
        ("terminal-nodes","place query to a terminal node")
        ("one-placement-only", "place only once despite equally good hits")
        ("454", "correct homopolymer error (DNA)")
        ("guided","guided placement with TID tags")
        ("pileup-alignment","make pileup alignment")

    ;

    boost::program_options::options_description reads_alignment2("Additional alignment extension options",100);
    reads_alignment2.add_options()
        ("homopolymer", "correct homopolymer error (more agressively)")
        ("pacbio","correct for missing data in PacBio reads (DNA)")
        ("query-distance", po::value<float>()->default_value(0.1,"0.1"), "evolutionary distance from pseudo-root")
        ("min-query-overlap", po::value<float>()->default_value(0.5,"0.5"), "overlap threshold for query and reference")
        ("overlap-with-any","accept query overlap with any sequence")
        ("min-query-identity", po::value<float>()->default_value(0.5,"0.5"), "identity threshold for aligned sites")
        ("pair-read-gap-extension", po::value<float>(), "paired read spacer extension probability (DNA)")
        ("score-only-ungapped","score query placement only on ungapped sites")
        ("score-ungapped-limit",po::value<float>()->default_value(0.1,"0.1"),"max. ungapped proportion")
        ("output-discarded-queries","output discarded queries to a file")
    ;

    boost::program_options::options_description reads_alignment3("Alignment extension output options",100);
    reads_alignment3.add_options()
        ("prune-extended-alignment","remove closely related sequences")
        ("prune-keep-number",po::value<int>()->default_value(0),"prune output and keep N most distantly related sequences")
        ("prune-keep-threshold",po::value<float>(),"prune output and remove sequences with distance below threshold")
        ("prune-keep-closest","prune output and keep only closest reference sequences")
        ("trim-extended-alignment","remove terminal reference sequences")
        ("trim-keep-sites",po::value<int>()->default_value(15),"trim distance around queries")
        ("use-consensus", "use consensus for query ancestors")
        ("show-contig-ancestor", "fill contig gaps with ancestral sequence")
        ("consensus-minimum", po::value<int>()->default_value(5), "threshold for inclusion in contig")
        ("consensus-minimum-proportion", po::value<float>()->default_value(0.5,"0.5"), "threshold for inclusion in contig")
    ;

    boost::program_options::options_description pileup("Pileup alignment options",100);
    pileup.add_options()
        ("min-orf-length", po::value<int>()->default_value(50),"minimum ORF length to be considered (DNA)")
        ("min-orf-coverage", po::value<float>(),"minimum ORF coverage to be considered (DNA)")
        ("query-cluster-attempts", po::value<int>()->default_value(1),"attempts to find overlap")
        ("build-contigs", "build contigs of query clusters")
        ("inlude-parent-in-contig", "include also ancestral parent in contigs")
        ("use-duplicate-weights", "use NumDuplicates=# to weight consensus counts")
        ("no-read-ordering","do not order reads by duplicate number")
        ("output-consensus", "output contig consensus alone")
    ;


    boost::program_options::options_description reads_quality("Quality options (FASTQ)",100);
    reads_quality.add_options()
        ("no-fastq", "do not use Q-scores")
        ("qscore-minimum", po::value<int>()->default_value(10), "threshold to mask low Q-score sites")
        ("perfect-reference", "assume perfect reference alignment")
    ;

    boost::program_options::options_description anchoring("Anchoring options (Exonerate)",100);
    anchoring.add_options()
        ("no-anchors","no anchoring (default: use Exonerate to anchor alignment)")
        ("exonerate-hit-length", po::value<int>()->default_value(30), "Exonerate hit length for anchor (default)")
        ("exonerate-hit-trim", po::value<int>()->default_value(5), "Exonerate hit trim length")
        ("exonerate-hit-score", po::value<int>(), "Exonerate hit score for anchor")
        ("anchors-offset", po::value<int>()->default_value(15), "anchors offset for alignment")
        ("anchoring-threshold",po::value<float>()->default_value(1,"1.0"),"anchoring coverage threshold for skipping")
        ("use-prefix-anchors","use prefix approach to anchor alignment")
        ("prefix-hit-length", po::value<int>()->default_value(30), "prefix hit length for anchor")
        ("keep-temp-files","keep temporary files (mainly for debugging)")
    ;

    boost::program_options::options_description exonerate("Fast placement options (Exonerate)",100);
    exonerate.add_options()
        ("exhaustive-placement","if Exonerate fails, use PAGAN to place the query")
        ("own-placement","use only PAGAN to place the query")
        ("use-exonerate-local","use Exonerate local to map queries to nodes")
        ("exonerate-local-keep-best",po::value<int>()->default_value(6),"keep best # of local matches")
        ("exonerate-local-keep-above",po::value<float>(),"keep local matches above #% of the best score")
        ("use-exonerate-gapped","use Exonerate gapped to map queries to nodes")
        ("exonerate-gapped-keep-best",po::value<int>()->default_value(3),"keep best # of gapped matches")
        ("exonerate-gapped-keep-above",po::value<float>(),"keep gapped matches above #% of the best score")
        ("keep-despite-exonerate-fails", "keep queries that Exonerate fails to align")
    ;

#ifdef NCBI_TOOLKIT
    boost::program_options::options_description ncbitoolkit("Anchoring options (NCBI toolkit)",100);
    ncbitoolkit.add_options()
        ("ncbi-threshold-overlap-total", po::value<int>()->default_value(50), "distance for fully overlapping hits to be accepted as anchors")
        ("ncbi-threshold-overlap-partly", po::value<int>()->default_value(400), "distance for partially overlapping hits to be accepted as anchors (should be higher than total overlapping)")
        ("blast-wordsize", po::value<int>()->default_value(-1), "wordsize for BLAST (default value: 11 for nucleic acids, 3 for amino acids)")
        ("blast-word-threshold", po::value<int>()->default_value(-1), "word threshold for BLAST (amino acids only, default value: 11)")
        ("blast-match-reward", po::value<int>()->default_value(-1), "match reward for BLAST (nucleic acids only, default value: 2)")
        ("blast-mismatch-penalty", po::value<int>()->default_value(999), "mismatch penalty for BLAST (nucleic acids only, default value: -3)")
        ("blast-scoring-matrix", po::value<string>()->default_value("BLOSUM62"), "scoring matrix for BLAST (amino acids only, possible values: BLOSUM45, BLOSUM50, BLOSUM62 (default), BLOSUM80, BLOSUM90)")
        ("memory-for-single-alignment", po::value<int>()->default_value(4000), "megabytes of memory allowed to use for a single alignment (using multiple threads may result mutiple concurrent alignment calculations)")
        ("force-gap", "force gaps in poorly prealigned areas to reduce memory consumption when exceeding memory limits")
        ("force-gap-threshold", po::value<int>()->default_value(40000), "min threshold size (=height*length | default=40000) for emply blocks in tunnel to be removed when exceeding memory limit")
        ("force-gap-wide-tunnel", "use wide tunnel when removing blocks (possibly results in fragmented gaps)")
    ;
#else
    boost::program_options::options_description ncbitoolkit("",100);
#endif


    boost::program_options::options_description obscure("Additional obscure options",100);
    obscure.add_options()
        ("rank-reads-for-nodes","rank reads within nodes for alignment")
        ("align-reads-at-root", "ignore tags and align reads at root")
        ("align-bad-reads-at-root", "align non-matching reads at root")
        ("use-identity-score", "choose target based on identity score")
        ("use-target-normalised-score", "choose target based on target-normalised substitution score")
        ("old-placement","old placement")
    ;

    boost::program_options::options_description graph("Graph options",100);
    graph.add_options()
        ("weight-sampled-edges", "use posterior scores to weight sampled edges")
        ("no-weight-transform", "no weight transform for sampled edges")
        ("cuberoot-weight-transform", "cuberoot weight transform for sampled edges")
    ;

    boost::program_options::options_description model("DNA/Protein model options",100);
    model.add_options()
        ("indel-rate", po::value<float>(), "insertion-deletion rate")
        ("gap-extension", po::value<float>(), "gap extension probability")
        ("end-gap-extension", po::value<float>(), "terminal gap extension probability")
        ("dna-kappa", po::value<float>(), "kappa")
        ("dna-rho", po::value<float>(), "rho")
        ("codons", "translate and align codons")
        ("use-aa-groups", "reconstruct amino-acid parsimony with 51 groups")
    ;

    boost::program_options::options_description tree_edit("Tree manipulation options",100);
    tree_edit.add_options()
        ("scale-branches", po::value<float>(), "scale tree branches")
        ("truncate-branches", po::value<float>()->default_value(0.2,"0.2"), "truncate tree branches")
        ("real-branches", "use real tree branch lengths")
        ("fixed-branches", po::value<float>(), "fixed length for tree branches")
        ("min-branch-length", po::value<float>(), "minimum length for tree branches")
    ;

    boost::program_options::options_description alignment("Alignment model options",100);
    alignment.add_options()
        ("any-skips-confirm-insertion", po::value<int>(), "#skips to confirm as insertion")
        ("match-skips-confirm-insertion", po::value<int>(), "#skips from match sites to confirm as insertion")
        ("branch-length-confirm-insertion", po::value<float>(), "total branch length skipped to confirm as insertion")
        ("branch-skip-weight-per-distance", po::value<float>(), "weighted (by branch length unit) probability for site(s) being skipped over and later matched (>default<)")
        ("branch-skip-penalty-per-branch", po::value<float>(), "fixed probability for site(s) being skipped over and later matched")
        ("keep-all-edges","nothing of those -- keep everything forever")
    ;

    boost::program_options::options_description output("Graph output options",100);
    output.add_options()
        ("mpost-graph-file", po::value<string>(), "sequence graph file for metapost")
        ("output-alignment-graphs", "include aligned graphs")
        ("output-leaf-graphs", "include terminal sequences")
        ("mpost-posterior-plot-file", po::value<string>(), "posterior plot file for metapost")
        ("plot-slope-up", "plot viterbi path climbing up")
    ;

    boost::program_options::options_description debug("Debugging and testing options",100);
    debug.add_options()
        ("output-nhx-tree", "output alignment tree (with NHX TID tags)")
        ("output-ancestors", "include ancestors in outfile")
        ("test-every-node","test every node for each query")
        ("test-every-internal-node","test every internal node for each query")
        ("test-every-terminal-node","test every terminal node for each query")
        ("no-preselection","do not preselect targets with Exonerate")
        ("score-as-dna", "score protein/ORFs as DNA (translated placement)")
        ("mostcommon", "use mostcommon for ambiguity")
        ("full-probability", "compute full probability")
        ("output-graph","output ancestral graph")
        ("sample-path", "sample the alignment path from posterior probabilities")
        ("ins-rate", po::value<float>(), "insertion rate (per substitution)")
        ("del-rate", po::value<float>(), "deletion rate (per substitution)")
        ("check-valid-graphs", "check that fwd and bwd edges are identical")
        ("full-help", "display full-help message")
        ("ambiguity-factor", po::value<float>(), "multiplier for subst. score of ambiguity characters")
        ("no-log-odds", "do not use log-odds substitutions scores")
        ("time", "track time (debugging)")
        ("recompute-reference-alignment-model", "recompute reference alignment model")
        ("no-score-scaling","no subsistitution score scaling")
        ("plot-anchors-for-R","plot for R")
        ("no-reduced-terminal-penalties", "no reduced terminal penalties")
        ("hmmer-anchors","hmmer anchors")
        ("tid-for-subroot","placement at subroot only (for assembly)")
        ("assembly","placement at subroot only (for assembly)")
        ("boost","multi-threading with boost")
        ("quick","quick reconstruction")
        ("docker","believe that binaries are there")
    ;

    boost::program_options::options_description broken("Broken options",100);
    broken.add_options()
        ("sample-additional-paths", po::value<int>()->default_value(0), "sample additional paths from posterior probabilities")
    ;

    po::positional_options_description pd;
    pd.add("config-file", 1);

    full_desc.add(minimal).add(help_update).add(generic).add(generic2).add(reads_alignment).add(reads_alignment2).add(reads_alignment3).add(pileup).add(reads_quality)
            .add(anchoring).add(ncbitoolkit).add(exonerate).add(obscure).add(graph).add(model).add(tree_edit).add(alignment).add(output).add(debug).add(broken);

    max_desc.add(minimal).add(generic).add(generic2).add(reads_alignment).add(reads_alignment2).add(reads_alignment3).add(pileup).add(reads_quality)
            .add(anchoring).add(ncbitoolkit).add(exonerate).add(obscure).add(model).add(graph).add(tree_edit).add(alignment).add(output).add(help_update);

    desc.add(minimal).add(generic).add(reads_alignment).add(reads_alignment2).add(reads_quality).add(anchoring).add(ncbitoolkit).add(exonerate).add(pileup).add(model).add(tree_edit).add(alignment).add(help_update);

    min_desc.add(minimal).add(generic).add(reads_alignment).add(help_update);


//    po::store(po::parse_command_line(argc, argv, full_desc), vm);
    po::store(po::command_line_parser(argc, argv).options(full_desc).positional(pd).run(), vm);

    if(Settings::is("config-file"))
    {
        stringstream ss;
        ss<<"\nReading command line options from file '" << Settings::get("config-file").as<string>()<<"'.\n";
        Log_output::write_out(ss.str(),0);

        ifstream cfg(Settings::get("config-file").as<string>().c_str());

        if (!cfg)
        {
            Settings::info_noexit();

            stringstream ss;
            ss<<"Config file '" << Settings::get("config-file").as<string>()<<"' not found. Exiting.\n\n";
            Log_output::write_out(ss.str(),0);

            exit(1);
        }


        po::store(po::parse_config_file(cfg, full_desc), vm);
    }


    po::notify(vm);


    ////////////////////////////////////////////////////////////////////////////////////////////

    if(is("noise"))
        noise = get("noise").as<int>();

    if(is("silent"))
        noise = -1;

    if(is("exonerate-local-keep-best") && get("exonerate-local-keep-best").as<int>() > 0 )
        exonerate_local_keep_best = get("exonerate-local-keep-best").as<int>();

    if(is("use-exonerate-local") && get("exonerate-local-keep-best").as<int>() > 0 )
        exonerate_local_keep_best = get("exonerate-local-keep-best").as<int>();

    if(is("exonerate-gapped-keep-best") && get("exonerate-gapped-keep-best").as<int>() > 0 )
        exonerate_gapped_keep_best = get("exonerate-gapped-keep-best").as<int>();

    if(is("use-exonerate-gapped") && get("exonerate-gapped-keep-best").as<int>() > 0 )
        exonerate_gapped_keep_best = get("exonerate-gapped-keep-best").as<int>();


    if(is("very-fast-placement"))
    {
        exonerate_local_keep_best = 1;
        exonerate_gapped_keep_best = 0;

        if( is("use-exonerate-local") || /*is("exonerate-local-keep-best") ||*/ is("exonerate-local-keep-above") ||
            is("use-exonerate-gapped") || /*is("exonerate-gapped-keep-best") ||*/ is("exonerate-gapped-keep-above") ||
            is("exhaustive-placement") || is("keep-despite-exonerate-fails") )
        {
            Log_output::write_out("Please disable Exonerate-related options if using '--very-fast-placement'. Exiting.\n",0);
            exit(1);
        }
    }

    if(is("fast-placement"))
    {
        exonerate_local_keep_best = 5;
        exonerate_gapped_keep_best = 1;

        if( is("use-exonerate-local") || is("exonerate-local-keep-above") ||
            is("use-exonerate-gapped") || is("exonerate-gapped-keep-above") ||
            is("exhaustive-placement") || is("keep-despite-exonerate-fails") )
        {
            Log_output::write_out("Please disable Exonerate-related options if using '--fast-placement'. Exiting.\n",0);
            exit(1);
        }
    }

    if( (is("very-fast-placement") || is("fast-placement")) &&
            !( is("test-every-node") || is("test-every-internal-node") || is("test-every-terminal-node") || is("all-nodes") || is("internal-nodes") || is("terminal-nodes") ) )
    {
        Log_output::write_out("\nWarning: When using option '--(very-)fast-placement', either option '--[all|internal|terminal]-nodes'\nor a reference tree with TID tags should be used. "
                              "If the latter is true, this warning can be ignored.\n",0);
    }

    // this heuristic only works for placement
    //
    tunneling_coverage = get("anchoring-threshold").as<float>();

    if(not is("queryfile"))
    {
        tunneling_coverage = 1;
    }

    placement_target_nodes = Settings::all_nodes;

    if(is("test-every-internal-node") || is("internal-nodes"))
    {
        placement_target_nodes = Settings::internal_nodes;
    }
    else if(is("test-every-terminal-node") || is("terminal-nodes"))
    {
        placement_target_nodes = Settings::terminal_nodes;
    }
    else if(is("test-every-node") || is("all-nodes"))
    {
        placement_target_nodes = Settings::all_nodes;
    }


    if(is("own-placement"))
    {
        exonerate_local_keep_best = 0;
        exonerate_gapped_keep_best = 0;
    }

    if(is("no-preselection") || is("guided"))
    {
        placement_preselection = false;
        placement_target_nodes = Settings::tid_nodes;
    }

    /////////////////////////////////////////////////////////////////////

    if (vm.count("help")) {
        this->help();
        return 1;
    }

    if (vm.count("full-help")) {
        this->help_all();
        return 1;
    }

    if (vm.count("version")) {
        this->check_version();
        return 1;
    }



    if(Settings::is("config-log-file"))
    {
        stringstream ss;
        ss << endl<< "Writing command line options to file '" << Settings::get("config-log-file").as<string>()<<"'."<<endl;
        Log_output::write_out(ss.str(),1);

        ofstream log_out(Settings::get("config-log-file").as<string>().c_str());
        time_t s_time;
        time( &s_time );
        log_out <<this->print_log_msg()<< "#\n# Analysis started: " << asctime( localtime( &s_time ) );
        log_out<<"# Command line arguments:"<<endl<<endl;

        po::parsed_options opts = parse_command_line(argc, argv, full_desc);

        typedef vector< po::basic_option<char> > vec_opt;

        for(vec_opt::iterator iter = opts.options.begin(); iter != opts.options.end(); ++iter)
        {
            po::basic_option<char>& option = *iter;

            if(option.string_key == "config-log-file" || option.string_key == "config-file")
                continue;

            stringstream ss_opt;
            typedef vector< basic_string<char> > vec_string;

            for(vec_string::iterator s_iter = option.value.begin(); s_iter != option.value.end(); ++s_iter)
            {
                    ss_opt << *s_iter;
            }
            if(ss_opt.str().length()>0)
                log_out<<option.string_key << " = " << ss_opt.str()<<endl;
            else
                log_out<<option.string_key << " = 1"<<endl;
        }

        if(Settings::is("config-file"))
        {
            log_out<<"\n# Additional arguments from file '"<<Settings::get("config-file").as<string>()<<"':"<<endl<<endl;

            ifstream cfg(Settings::get("config-file").as<string>().c_str());
            opts = po::parse_config_file(cfg, full_desc);

            for(vec_opt::iterator iter = opts.options.begin(); iter != opts.options.end(); ++iter)
            {
                po::basic_option<char>& option = *iter;

                if(option.string_key == "config-log-file" || option.string_key == "config-file")
                    continue;

                stringstream ss_opt;
                typedef vector< basic_string<char> > vec_string;

                for(vec_string::iterator s_iter = option.value.begin(); s_iter != option.value.end(); ++s_iter)
                {
                        ss_opt << *s_iter;
                }
                if(ss_opt.str().length()>0)
                    log_out<<option.string_key << " = " << ss_opt.str()<<endl;
                else
                    log_out<<option.string_key << " = 1"<<endl;
            }
        }
        log_out<<endl;
    }

    return 0;
}

void Settings::print_msg()
{
    stringstream ss;
    ss<<"\nPAGAN2 v."<<version<<" ("<<date<<"). (C) 2010-2019 by Ari Löytynoja <ari.loytynoja@gmail.com>.\n";
    ss<<" This program is provided \"as-is\", with NO WARRANTY whatsoever; this is a development version\n and may contain bugs.\n";
    Log_output::write_out(ss.str(),0);
}

string Settings::print_log_msg()
{
    stringstream tmp;
    tmp<<"\n# PAGAN2 v."<<version<<" ("<<date<<"). (C) 2010-2019 by Ari Löytynoja <ari.loytynoja@gmail.com>.\n";
    tmp<<"# This program is provided \"as-is\", with NO WARRANTY whatsoever; this is a development version\n# and may contain bugs.\n";
    return tmp.str();
}

void Settings::help()
{
    this->print_msg();
    stringstream ss;
    ss<< desc << "\n";
    Log_output::write_out(ss.str(),0);
    exit(0);
}

void Settings::help_all()
{
    this->print_msg();
    stringstream ss;
    ss<< max_desc << "\n";
    Log_output::write_out(ss.str(),0);
    exit(0);
}

void Settings::info()
{
    this->print_msg();
    stringstream ss;
    ss<< min_desc << "\n";
    Log_output::write_out(ss.str(),0);
    exit(0);
}

void Settings::info_noexit()
{
    this->print_msg();
    stringstream ss;
    ss<< min_desc << "\n\n";
    Log_output::write_out(ss.str(),0);
}

void Settings::check_version()
{
    Check_version cv(version);
    exit(0);
}

int   Settings::noise             = 0;
float Settings::resize_factor     = 1.5;

int   Settings::exonerate_local_keep_best = 0;
int   Settings::exonerate_gapped_keep_best = 0;

float Settings::tunneling_coverage = 1;

int Settings::placement_target_nodes = Settings::tid_nodes;

bool Settings::placement_preselection = true;
