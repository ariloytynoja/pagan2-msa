/***************************************************************************
 *   Copyright (C) 2010-2019 by Jouko Niinimäki & Ari Loytynoja            *
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

#ifndef NCBI_BLAST_H
#define NCBI_BLAST_H

#include <objects/seqalign/Dense_diag.hpp>
#include <objects/seqalign/Dense_seg.hpp>
#include <objects/seqloc/Seq_loc.hpp>
#include <objects/seqalign/Score.hpp>

#include <objmgr/object_manager.hpp>

#include <algo/blast/api/bl2seq.hpp>
#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/blast_options_handle.hpp>

#include <objmgr/scope.hpp>

#include <limits>
#include <sstream>
#include <queue>

#include "utils/substring_hit.h"
#include "utils/log_output.h"

namespace ppa {

USING_NCBI_SCOPE;
USING_SCOPE(objects);
using namespace ppa;
using namespace std;

enum mol_type{
    dna,
    rna,
    aminoa,
    nucla
};

enum strand{
    strand_plus,
    strand_minus,
    strand_both
};

enum score_matrix{
    blosum45,
    blosum50,
    blosum62,
    blosum80,
    blosum90
};

struct Alignment_set{
    string seq1_id;
    string seq2_id;
    vector<Substring_hit> hits;
    int total_score;
    bool score_strand_plus;

    Alignment_set() : total_score(0), score_strand_plus(false){

    }
    void print(stringstream& s){
        printSimple(s);

        for (vector<Substring_hit>::iterator hit = hits.begin(); hit != hits.end(); ++hit){
            s << "    ";
            hit->print(s);
        }
    }
    void printSimple(stringstream& s){
        s << "Alignment: " << seq1_id << " & " << seq2_id << ", hits: " << hits.size() << ", tot_score: " << total_score << (score_strand_plus ? "" : " (reverse)");
    }
    bool operator > (const Alignment_set& bh) const{
        return (total_score > bh.total_score);
    }
    bool operator < (const Alignment_set& bh) const{
        return (total_score < bh.total_score);
    }
    void countTotalScore(){
        int plus_score = 0;
        int minus_score = 0;

        for(vector<Substring_hit>::iterator hit = hits.begin(); hit != hits.end(); hit++){
            if(hit->plus_strand_1 && hit->plus_strand_2){
                plus_score += hit->score;
            }else{
                minus_score += hit->score;
            }
        }

        if(minus_score > plus_score){
            stringstream s;
            s << "Note: Align '" << seq1_id << " & " << seq2_id << "' has higher reverse strand score!";
            Log_output::write_warning(s.str(), 3);
            score_strand_plus = false;
        }else{
            score_strand_plus = true;
        }

        total_score = max(plus_score, minus_score);
    }
};

struct compare_alignment_sets{
    bool operator() ( const Alignment_set& first, const Alignment_set& second ) const{
        return first.total_score > second.total_score;
    }
};
struct compare_alignment_set_ptrs{
    bool operator() ( const Alignment_set* first, const Alignment_set* second ) const{
        return first->total_score > second->total_score;
    }
};

struct Blast_options{
    int wordsize;
    int word_threshold_score;
    int match_reward;
    int mismatch_penalty;
    score_matrix scoring_matrix;
    strand strand_opt;


    Blast_options() : wordsize(-1), word_threshold_score(-1), match_reward(-1), mismatch_penalty(999), scoring_matrix(blosum62), strand_opt(strand_plus){

    }
    bool operator == (const Blast_options& o){
        if(wordsize == o.wordsize && word_threshold_score == o.word_threshold_score && match_reward == o.match_reward && mismatch_penalty == o.mismatch_penalty && scoring_matrix == o.scoring_matrix && strand_opt == o.strand_opt){
            return true;
        }
        return false;
    }
    bool operator != (const Blast_options& o){
        return !(*this == o);
    }
    bool isDefault(){
        if(wordsize == -1 && word_threshold_score == -1 && match_reward == -1 && mismatch_penalty == 999 && scoring_matrix == blosum62 && strand_opt == strand_plus){
            return true;
        }
        return false;
    }
};

struct Alignment_bank{
    map<string, vector<Alignment_set> > alignment_sets;
    map<string, string> querys;
    vector<string> target_ids;

    mol_type mol;
    Blast_options options;
};

typedef std::priority_queue< Alignment_set*,  std::vector<Alignment_set*>, compare_alignment_set_ptrs > Alignment_set_queue;


class NCBIBlaster{

private:
    void setDefaultOptions(CRef<blast::CBlastOptionsHandle>& options, mol_type molecyle);
    blast::SSeqLoc* getSeqFromString(const string& id, const string& sequence, mol_type molecyle);
    blast::SSeqLoc* getSeqLocFromSeqEntry(CRef<CSeq_entry>& entry);
    void setOptions(CRef<blast::CBlastOptionsHandle>& options_handle, Blast_options& options, mol_type molecyle);
    void parseSingleResults(blast::TSeqAlignVector& results, vector<Substring_hit>& hits);
    void parseMultiResults(blast::TSeqAlignVector& results, vector<Alignment_set>& align_sets);
    void parseDiagAlignment(CRef<CSeq_align>& seq_align, vector<Substring_hit>& hits);

    int executeBlast(blast::TSeqLocVector& querys, blast::TSeqLocVector& subjects, blast::TSeqAlignVector& alignments, mol_type mol, Blast_options* options = 0);

public:
    virtual int runSingleBlast(string& query, string& subject, string& query_id, string& subject_id, vector<Substring_hit>& hits, mol_type mol, Blast_options* options = 0);
    int runMultiBlast(map<string, string>& querys, map<string, string>& subjects, vector<Alignment_set>& align_sets, mol_type mol, Blast_options* options = 0);
    NCBIBlaster();


};

class BankBlaster : NCBIBlaster{
private:

    static Alignment_bank bank;

    bool sameOptions(mol_type mol, Blast_options* options){
        if(mol == bank.mol){
            if(options == 0 && bank.options.isDefault()){
                return true;
            }
            if(options != 0 && *options == bank.options){
                return true;
            }
        }
        return false;
    }
    void setBankBlastOptions(mol_type mol, Blast_options* options = 0){
        bank.mol = mol;
        if(options != 0){
            bank.options = *options;
        }
    }
    void printBank(stringstream& s){
        for(map<string, string>::iterator query = bank.querys.begin(); query != bank.querys.end(); query++){
            s << "  query: " << query->first << endl;
            vector<Alignment_set>& sets_for_q = bank.alignment_sets[query->first];

            for(vector<Alignment_set>::iterator aset = sets_for_q.begin(); aset != sets_for_q.end(); aset++){
                s << "    ";
                aset->printSimple(s);
                s << endl;
            }
        }
    }

public:
    int runSingleBlast(string& query, string& subject, string& query_id, string& subject_id, vector<Substring_hit>& hits, mol_type mol, Blast_options* options = 0);
    void initBankAlignments(map<string, string>& querys, map<string, string>& targets, int save_x_best, mol_type mol, Blast_options* options = 0);
    void insertOneTarget(const string& id, const string& sequence);
    void insertTargets(map<string, string>& targets);
    void pickXBestAlignsFromBank(vector<Alignment_set*>& best_aligns, string query_id, int x);
    void copyXBestAlignsFromBank(vector<Alignment_set>& best_aligns, string query_id, int x);
    bool isBestAlignReverse(string& query_id);

};



}

#endif
