#include "ncbi_blast.h"

USING_NCBI_SCOPE;
USING_SCOPE(objects);
using namespace ppa;

Alignment_bank BankBlaster::bank;

/**
 * @brief NCBIBlaster::NCBIBlaster
 */
NCBIBlaster::NCBIBlaster(){

}

/**
 * @brief NCBIBlaster::setDefaultOptions constructs default options for blast based on given molecyle
 * @param options
 * @param molecyle
 */
void NCBIBlaster::setDefaultOptions(CRef<blast::CBlastOptionsHandle>& options, mol_type molecyle){
    if(molecyle == aminoa){
        options.Reset(blast::CBlastOptionsFactory::Create(ncbi::blast::eBlastp));
        options->SetOptions().SetCompositionBasedStats(eNoCompositionBasedStats);
    }else{
        options.Reset(blast::CBlastOptionsFactory::Create(ncbi::blast::eBlastn));
    }

    options->SetGappedMode(false);
    options->SetOptions().SetStrandOption(eNa_strand_plus);
}

void NCBIBlaster::setOptions(CRef<blast::CBlastOptionsHandle>& options_handle, Blast_options& options, mol_type molecyle){

    if(options.wordsize >= 0){
        options_handle->SetOptions().SetWordSize(options.wordsize);
    }
    if(options.word_threshold_score >= 0){
        options_handle->SetOptions().SetWordThreshold(options.word_threshold_score);
    }

    switch(options.strand_opt){
        case strand_plus:
            options_handle->SetOptions().SetStrandOption(eNa_strand_plus);
            break;
        case strand_minus:
            options_handle->SetOptions().SetStrandOption(eNa_strand_minus);
            break;
        case strand_both:
            options_handle->SetOptions().SetStrandOption(eNa_strand_both);
            break;
    }

    if(molecyle == aminoa){
        switch(options.scoring_matrix){
            case blosum45:
                options_handle->SetOptions().SetMatrixName("BLOSUM45");
                break;
            case blosum50:
                options_handle->SetOptions().SetMatrixName("BLOSUM50");
                break;
            case blosum62:
                options_handle->SetOptions().SetMatrixName("BLOSUM62");
                break;
            case blosum80:
                options_handle->SetOptions().SetMatrixName("BLOSUM80");
                break;
            case blosum90:
                options_handle->SetOptions().SetMatrixName("BLOSUM90");
                break;
        }
    }else{
        if(options.match_reward >= 0){
            options_handle->SetOptions().SetMatchReward(options.match_reward);
        }
        if(options.mismatch_penalty < 999){
            options_handle->SetOptions().SetMismatchPenalty(options.mismatch_penalty);
        }
    }

}

/**
 *
 * @brief NCBIBlaster::getSeqFromString returns new SSeqLoc* for given sequence string
 * @param molecyle
 * @param id
 * @param sequence
 * @return
 */
blast::SSeqLoc* NCBIBlaster::getSeqFromString(const string& id, const string& sequence, mol_type molecyle){

    CRef<CSeq_entry> seq_entry(new CSeq_entry);
    CRef<CBioseq> bioseq(new CBioseq);

    bioseq->SetInst().SetRepr(CSeq_inst::eRepr_raw);


    switch(molecyle){
        case dna:
            bioseq->SetInst().SetMol(CSeq_inst::eMol_dna);
            bioseq->SetInst().SetSeq_data().SetIupacna().Set(sequence);
            //bioseq->SetInst().SetStrand(CSeq_inst::eStrand_ds);
            break;
        case rna:
            bioseq->SetInst().SetMol(CSeq_inst::eMol_rna);
            bioseq->SetInst().SetSeq_data().SetIupacna().Set(sequence);
            break;
        case aminoa:
            bioseq->SetInst().SetMol(CSeq_inst::eMol_aa);
            bioseq->SetInst().SetSeq_data().SetIupacaa().Set(sequence);
            //bioseq->SetInst().SetStrand(CSeq_inst::eStrand_ss);
            break;
        case nucla: //nucleotic acid
            bioseq->SetInst().SetMol(CSeq_inst::eMol_na);
            bioseq->SetInst().SetSeq_data().SetIupacna().Set(sequence);
            break;
    }

    bioseq->SetInst().SetLength(sequence.length());


    bioseq->SetId().push_back(CRef<CSeq_id>(new CSeq_id(CSeq_id::e_Local, id, kEmptyStr)));

    seq_entry->SetSeq(*bioseq);

    return getSeqLocFromSeqEntry(seq_entry);
}
/**
 * @brief NCBIBlaster::getSeqLocFromSeqEntry creates and returns new SSeqLoc* for given CSeq_entry
 * @param entry
 * @return
 */
blast::SSeqLoc* NCBIBlaster::getSeqLocFromSeqEntry(CRef<CSeq_entry>& entry){

    const CSeq_id* seqid = entry->GetSeq().GetFirstId();

    //cout << "id: " << seqid->AsFastaString() << endl;

    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    CRef<CScope> scope(new CScope(*objmgr));

    scope->AddDefaults();
    scope->AddTopLevelSeqEntry(*entry);

    CBioseq_Handle seq_handle = scope->GetBioseqHandle(*seqid);

    CRef<CSeq_loc> seqloc(new CSeq_loc());

    seqloc->SetInt().SetTo(seq_handle.GetInst_Length() - 1);
    seqloc->SetInt().SetFrom(0);
    seqloc->SetInt().SetId().Assign(*entry->GetSeq().GetId().front());

    return new blast::SSeqLoc(seqloc, scope);
}
int NCBIBlaster::runSingleBlast(string& query, string& subject, string& query_id, string& subject_id, vector<Substring_hit>& hits, mol_type mol, Blast_options* options /*=0*/){

    int r = 0;
    #pragma omp critical(blast)
    {


        auto_ptr<blast::SSeqLoc> q_loc(getSeqFromString(query_id, query, mol));
        auto_ptr<blast::SSeqLoc> s_loc(getSeqFromString(subject_id, subject, mol));

        blast::TSeqLocVector q_vec;
        blast::TSeqLocVector s_vec;

        q_vec.push_back(*q_loc);
        s_vec.push_back(*s_loc);

        //vector<Alignment_set> align_sets;
        blast::TSeqAlignVector alignments;
        r = executeBlast(q_vec, s_vec, alignments, mol, options);
        parseSingleResults(alignments, hits);

        //hits = align_sets[0].hits;
    }
    return r;
}

int NCBIBlaster::runMultiBlast(map<string, string>& querys, map<string, string>& subjects, vector<Alignment_set>& align_sets, mol_type mol, Blast_options* options /*=0*/){

    int r = 0;
    #pragma omp critical(blast)
    {

        blast::TSeqLocVector q_vec;
        blast::TSeqLocVector s_vec;

        for(map<string, string>::iterator seq = querys.begin(); seq != querys.end(); ++seq){
            auto_ptr<blast::SSeqLoc> q_loc(getSeqFromString(seq->first, seq->second, mol));
            q_vec.push_back(*q_loc);
        }
        for(map<string, string>::iterator seq = subjects.begin(); seq != subjects.end(); ++seq){
            auto_ptr<blast::SSeqLoc> s_loc(getSeqFromString(seq->first, seq->second, mol));
            s_vec.push_back(*s_loc);
        }

        blast::TSeqAlignVector alignments;
        r = executeBlast(q_vec, s_vec, alignments, mol, options);
        parseMultiResults(alignments, align_sets);

    }
    return r;
}

/**
 * @brief NCBIBlaster::executeBlast runs blast with with two given sequences
 * @param query
 * @param subject
 * @param hits
 * @param mol
 * @return
 */
int NCBIBlaster::executeBlast(blast::TSeqLocVector& querys, blast::TSeqLocVector& subjects, blast::TSeqAlignVector& alignments, mol_type mol, Blast_options* options /*=0*/){

    CRef<blast::CBlastOptionsHandle> options_handle;
    setDefaultOptions(options_handle, mol);
    if(options != 0){
        setOptions(options_handle, *options, mol);
    }

    int r = 0;
    try {
        blast::CBl2Seq blaster(querys, subjects, *options_handle);
        //#pragma omp critical
        alignments = blaster.Run();

        //for (int i = 0; i < (int) alignments.size(); ++i)
          //  std::cout << MSerial_AsnText << *alignments[i];

    } catch (const CException& e) {
        cout << e.GetMsg() << endl;
        r = 1;
    }

    return r;
}
/**
 * @brief NCBIBlaster::parseResults
 * @param results
 * @param hits
 */
void NCBIBlaster::parseSingleResults(blast::TSeqAlignVector& results, vector<Substring_hit>& hits){

    //cout << "Tseqalignvectorin pituus: " << results.size() << endl;
    CRef<CSeq_align_set> resaligns = results[0];


    //cout << "Seq_align_setin pituus: " << resaligns->Get().size() << endl;
    if(resaligns.GetObject().CanGet()){
        ///list of alignments of one sequence pair
        list<CRef<CSeq_align> > l = resaligns->Get();
        for(list<CRef<CSeq_align> >::iterator seq_align_it = l.begin(); seq_align_it != l.end(); seq_align_it++){
            if((*seq_align_it)->GetType() == CSeq_align_Base::eType_diags){
                parseDiagAlignment(*seq_align_it, hits);
            }else{
                Log_output::write_warning("Blast result parse error: unsupported alignment type (this should not happen..). Alignment discarded.", 0);
            }//TODO: other types (especially partial type)
        }
    }

}

void NCBIBlaster::parseMultiResults(blast::TSeqAlignVector& results, vector<Alignment_set>& align_sets){

    //cout << "Tseqalignvectorin pituus: " << results.size() << endl;
    for(int i = 0; i < (int)results.size(); i++){
        CRef<CSeq_align_set> resaligns = results[i];

        Alignment_set align;
        align_sets.push_back(align);

        //cout << "Seq_align_setin pituus: " << resaligns->Get().size() << endl;
        if(resaligns.GetObject().CanGet()){
            ///list of alignments of one sequence pair
            list<CRef<CSeq_align> > l = resaligns->Get();
            for(list<CRef<CSeq_align> >::iterator seq_align_it = l.begin(); seq_align_it != l.end(); seq_align_it++){

                Alignment_set& current_align = align_sets.back();
                current_align.seq1_id = (*seq_align_it)->GetSeq_id(0).GetSeqIdString();
                current_align.seq2_id = (*seq_align_it)->GetSeq_id(1).GetSeqIdString();

                if((*seq_align_it)->GetType() == CSeq_align_Base::eType_diags){
                    parseDiagAlignment(*seq_align_it, current_align.hits);
                }else{
                    Log_output::write_warning("Blast result parse error: unsupported alignment type (this should not happen..). Alignment discarded.", 0);
                }//TODO: other types (especially partial type)
            }
        }
    }

}

/**
 * @brief NCBIBlaster::parseDiagAlignment fetches the relevant data from a diagonal alignment
 * @param hits
 * @param seq_align
 */
void NCBIBlaster::parseDiagAlignment(CRef<CSeq_align>& seq_align, vector<Substring_hit>& hits){

    ///segments object (that has list of segments)
    const CSeq_align_Base::C_Segs& segs = seq_align->GetSegs();
    if(segs.IsDendiag()){

        ///iterating through list of segments
        for(list<CRef<CDense_diag> >::const_iterator ddiag_it = segs.GetDendiag().begin(); ddiag_it != segs.GetDendiag().end(); ddiag_it++){

            Substring_hit hit;
            hit.score = 1;

            ///starting positions
            vector<unsigned int> starts = (*ddiag_it)->GetStarts();
            if(starts.size() == 2){
                hit.start_site_1 = (int) starts[0];
                hit.start_site_2 = (int) starts[1];
            }else{
                Log_output::write_warning("Blast result parse error: pairwise alignment (hit) does not have two start sites. Hit discarded.", 0);
                return;
            }

            ///strands
            vector<ENa_strand> strands = (*ddiag_it)->GetStrands();
            if(strands.size() == 2){
                switch(strands[0]){
                    case eNa_strand_plus:
                        hit.plus_strand_1 = true;
                        break;
                    case eNa_strand_minus:
                        hit.plus_strand_1 = false;
                        break;
                    case eNa_strand_unknown:
                        hit.plus_strand_1 = true;
                        break;
                    default:
                        Log_output::write_warning("Blast result parse error: invalid strand information. Hit discarded.", 0);
                        return;
                }
                switch(strands[1]){
                    case eNa_strand_plus:
                        hit.plus_strand_2 = true;
                        break;
                    case eNa_strand_minus:
                        hit.plus_strand_2 = false;
                        break;
                    case eNa_strand_unknown:
                        hit.plus_strand_2 = true;
                        break;
                    default:
                        Log_output::write_warning("Blast result parse error: invalid strand information. Hit discarded.", 0);
                        return;
                }
            }else{
                Log_output::write_warning("Blast result parse error: pairwise alignment (hit) has more or less than two strands. Hit discarded.", 0);
                return;
            }


            ///length
            hit.length = (*ddiag_it)->GetLen();

            ///scores
            vector<CRef<CScore> > scores = (*ddiag_it)->GetScores();
            for(int i = 0; i < (int) scores.size(); i++){
                const CObject_id& score_id = scores[i]->GetId();
                if(score_id.GetStr() == "score"){
                    hit.score = scores[i]->GetValue().GetInt();
                    break;

                }
            }


            hits.push_back(hit);
        }
    }
}





/**********************************************************************************************************************************************************/
/*                                                                    BankBlaster                                                                         */
/**********************************************************************************************************************************************************/

void printSimpleHits(vector<Substring_hit>& hits, string qid, string tid){
    Alignment_set a;
    a.hits = hits;
    a.seq1_id = qid;
    a.seq2_id = tid;
    a.countTotalScore();
    stringstream s;
    a.printSimple(s);
    Log_output::write_msg(s.str(), 0);
}

int BankBlaster::runSingleBlast(string& query, string& subject, string& query_id, string& subject_id, vector<Substring_hit>& hits, mol_type mol, Blast_options* options /*= 0*/){

    if(sameOptions(mol, options)){
        if(bank.alignment_sets.find(query_id) != bank.alignment_sets.end()){
            vector<Alignment_set>& asets = bank.alignment_sets[query_id];
            for(vector<Alignment_set>::iterator aset = asets.begin(); aset != asets.end(); aset++){
                if(aset->seq2_id == subject_id){
                    vector<Substring_hit>& aset_hits = aset->hits;
                    for(vector<Substring_hit>::iterator hit = aset_hits.begin(); hit != aset_hits.end(); hit++){
                        hits.push_back(*hit);
                    }
                    stringstream s;
                    #pragma omp critical(other)
                    {
                        s << "Alignment '" << query_id << " vs " << subject_id << "' found in bank: ";
                        aset->printSimple(s);
                        Log_output::write_msg(s.str(), 2);
                    }
                    return 0;
                }
            }
        }
    }

    stringstream s;
    #pragma omp critical(other)
    {
        s << "Running blast: '" << query_id << " vs " << subject_id << "' (not found in alignment bank)..";
        Log_output::write_msg(s.str(), 2);
    }
    int r = NCBIBlaster::runSingleBlast(query, subject, query_id, subject_id, hits, mol, options);

    return r;
}


/**
 * @brief BankBlaster::initBankAlignments
 * @param querys
 * @param targets
 * @param save_x_best
 */
void BankBlaster::initBankAlignments(map<string, string>& querys, map<string, string>& targets, int save_x_best, mol_type mol, Blast_options* options /*= 0*/){
    setBankBlastOptions(mol, options);
    vector<Alignment_set> alignments;
    runMultiBlast(querys, targets, alignments, bank.mol, &bank.options);

    map<string, vector<Alignment_set*> > order_sets;

    bank.querys = querys;
    for(map<string, string>::iterator q = querys.begin(); q != querys.end(); q++){
        //this->bank.querys.push_back(q->first);
        bank.alignment_sets.insert(pair<string, vector<Alignment_set> >(q->first, vector<Alignment_set>()));
        order_sets.insert(pair<string, vector<Alignment_set*> >(q->first, vector<Alignment_set*>()));
    }
    for(map<string, string>::iterator t = targets.begin(); t != targets.end(); t++){
        bank.target_ids.push_back(t->first);
    }

    ///copies alignment_set ptrs to map for ordering
    for(vector<Alignment_set>::iterator as = alignments.begin(); as != alignments.end(); as++){

        as->countTotalScore();

        if(as->total_score > 0){
            order_sets[as->seq1_id].push_back(&(*as));
        }
    }

    ///orders map and copies best sets to bank
    for(map<string, vector<Alignment_set*> >::iterator sets_for_q = order_sets.begin(); sets_for_q != order_sets.end(); sets_for_q++){
        string query_id = sets_for_q->first;
        vector<Alignment_set*>& asets = sets_for_q->second;

        sort(asets.begin(), asets.end(), compare_alignment_set_ptrs());

        vector<Alignment_set>& asets_for_q_final = bank.alignment_sets[query_id];

        vector<Alignment_set*>::iterator aset = asets.begin();
        int last_score = (*aset)->total_score;
        int index = save_x_best;
        while(aset != asets.end() && (last_score == (*aset)->total_score || index > 0)){

            asets_for_q_final.push_back(*(*aset));
            last_score = (*aset)->total_score;
            aset++;
            index--;
        }

        /*
        for(int i = 0; i < asets.size(); i++){
            asets_for_q_final.push_back(*asets[i]);
        }*/
    }

    stringstream s;
    s << "Alignment bank initiated with following alignment sets: " << endl;
    printBank(s);
    Log_output::write_msg(s.str(), 2);

}

void BankBlaster::insertOneTarget(const string& id, const string& sequence){
    map<string, string> t;
    t.insert(pair<string, string>(id, sequence));
    insertTargets(t);
}


/**
 * @brief BankBlaster::insertTargets
 * @param targets
 */
void BankBlaster::insertTargets(map<string, string>& targets){
    vector<Alignment_set> alignments;
    runMultiBlast(bank.querys, targets, alignments, bank.mol, &bank.options);

    ///add target ids to target vector
    for(map<string, string>::iterator t = targets.begin(); t != targets.end(); t++){
        bank.target_ids.push_back(t->first);
    }

    ///copies alignment_sets to map
    for(vector<Alignment_set>::iterator as = alignments.begin(); as != alignments.end(); as++){

        as->countTotalScore();

        if(as->total_score > 0){
            bank.alignment_sets[as->seq1_id].push_back((*as));
        }
    }
    ///reorders alignment_sets for each query
    for(map<string, string>::iterator query = bank.querys.begin(); query != bank.querys.end(); query++){
        vector<Alignment_set>& sets_for_q = bank.alignment_sets[query->first];
        sort(sets_for_q.begin(), sets_for_q.end(), compare_alignment_sets());
    }

    stringstream s;
    s << "Alignment bank updated and contains now following alignment sets: " << endl;
    printBank(s);
    Log_output::write_msg(s.str(), 3);
}

/**
 * @brief BankBlaster::pickXBestAlignsFromBank returns a vector of pointers to x best aligns for given query in bank
 * @param best_aligns
 * @param query_id
 * @param x
 */
void BankBlaster::pickXBestAlignsFromBank(vector<Alignment_set*>& best_aligns, string query_id, int x){

    if (bank.alignment_sets.find(query_id) == bank.alignment_sets.end()){
        stringstream s;
        s << "Warning: Query '" << query_id << "'' not found from alignment bank. No aligns seleceted.";
        Log_output::write_warning(s.str(), 0);
        return;
    }
    vector<Alignment_set>& asets_for_q = bank.alignment_sets[query_id];

    int i = 0;
    int last_score = 0;
    while(i < (int)asets_for_q.size() && ( i < x || last_score == asets_for_q[i].total_score)){
        best_aligns.push_back(&asets_for_q[i]);
        last_score = asets_for_q[i].total_score;
        i++;
    }


    for(int i = 0; i < x && i < (int)asets_for_q.size(); i++){

    }
}

/**
 * @brief BankBlaster::copyXBestAlignsFromBank creates a copies of x best aligns for given query in bank
 * @param best_aligns
 * @param query_id
 * @param x
 */
void BankBlaster::copyXBestAlignsFromBank(vector<Alignment_set>& best_aligns, string query_id, int x){

    if (bank.alignment_sets.find(query_id) == bank.alignment_sets.end()){
        stringstream s;
        s << "Query '" << query_id << "'' not found from alignment bank. No aligns seleceted.";
        Log_output::write_warning(s.str(), 0);
        return;
    }
    vector<Alignment_set>& asets_for_q = bank.alignment_sets[query_id];

    int i = 0;
    int last_score = 0;
    while(i < (int)asets_for_q.size() && ( i < x || last_score == asets_for_q[i].total_score)){
        best_aligns.push_back(asets_for_q[i]);
        last_score = asets_for_q[i].total_score;
        i++;
    }


    for(int i = 0; i < x && i < (int)asets_for_q.size(); i++){

    }
}

/**
 * @brief BankBlaster::isBestAlignReverse
 * @param query_id
 * @return
 */
bool BankBlaster::isBestAlignReverse(string& query_id){

    if(bank.alignment_sets.find(query_id) != bank.alignment_sets.end()){
        vector<Alignment_set>& asets = bank.alignment_sets[query_id];
        if(asets.size() > 0 && !asets[0].score_strand_plus){
            return true;
        }
    }
    return false;
}

