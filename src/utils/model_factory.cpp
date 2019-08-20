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

#include <iomanip>
#include <iostream>
#include <vector>
#include <ctime>
#include "utils/model_factory.h"
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/eigen.h"
#include "utils/log_output.h"

using namespace ppa;
using namespace std;

Model_factory::Model_factory(int s)
{
    sequence_data_type = s;

    charPi=0;
    charU=0;
    charV=0;
    charRoot=0;

    parsimony_table=0;
    child_parsimony_table=0;
    mostcommon_table = 0;
    char_ambiguity=0;

//    prev_distance = -1;
//    prev_is_local_alignment = false;

    if(sequence_data_type == Model_factory::dna)
        if(Settings_handle::st.is("codons"))
            this->define_codon_alphabet();
        else
            this->define_dna_alphabet();
    else if(sequence_data_type == Model_factory::protein)
        if(Settings_handle::st.is("use-aa-groups"))
            this->define_protein_alphabet_groups();
        else
            this->define_protein_alphabet();
    else if(sequence_data_type == Model_factory::codon)     
        this->define_codon_alphabet();
    else
    {
        Log_output::write_out("Model_factory(): invalid sequence data type. Exiting\n",0);
        exit(1);
    }
}

Model_factory::~Model_factory()
{
    if(charPi!=0)
        delete charPi;

    if(charU!=0)
        delete charU;

    if(charV!=0)
        delete charV;

    if(charRoot!=0)
        delete charRoot;

    if(parsimony_table!=0)
        delete parsimony_table;

    if(child_parsimony_table!=0)
        delete child_parsimony_table;

    if(mostcommon_table!=0)
        delete mostcommon_table;

    if(char_ambiguity!=0)
        delete char_ambiguity;

//    if(model!=0)
//        delete model;
}

/*******************************************/

string Model_factory::dna_char_alphabet = "ACGT";
string Model_factory::dna_full_char_alphabet = "ACGTRYMKWSBDHVN";
string Model_factory::protein_char_alphabet = "ARNDCQEGHILKMFPSTWYV";
string Model_factory::protein_full_char_alphabet = "";

vector<string> Model_factory::dna_character_alphabet = vector<string>();
vector<string> Model_factory::dna_full_character_alphabet = vector<string>();
vector<string> Model_factory::protein_character_alphabet = vector<string>();
vector<string> Model_factory::protein_full_character_alphabet = vector<string>();
vector<string> Model_factory::codon_character_alphabet = vector<string>();
vector<string> Model_factory::codon_full_character_alphabet = vector<string>();

vector<string> Model_factory::ancestral_character_alphabet= vector<string>();


/*
 * Definition of DNA alpahbet including ambiguity characters.
 */
void Model_factory::define_dna_alphabet()
{
    char_alphabet = this->dna_char_alphabet;
    full_char_alphabet = this->dna_full_char_alphabet;

    char_as = char_alphabet.length();
    char_fas = full_char_alphabet.length();

    int n_residues[] = {1,1,1,1,2,2,2,2,2,2,3,3,3,3,4};
    string ambiguity[] = {"A","C","G","T","AG","CT","AC","GT","AT","CG","CGT","AGT","ACT","ACG","ACGT"};

    for(int i=0;i<15;i++)
    {
        Char_symbol letter;

        letter.index = i;
        letter.symbol = full_char_alphabet.at(i);
        letter.n_units = n_residues[i];
        for(int j=0;j<letter.n_units;j++)
            letter.residues.push_back(ambiguity[i].at(j));
        char_symbols.push_back(letter);
    }

    /*
     * Table for resolving the parental state from the parsimony ambiguity alphabet.
     */

    int a = 1;
    int c = 2;
    int g = 4;
    int t = 8;

    int r = (a|g);
    int y = (c|t);
    int m = (a|c);
    int k = (g|t);
    int w = (a|t);
    int s = (c|g);

    int b = (c|g|t);
    int d = (a|g|t);
    int h = (a|c|t);
    int v = (a|c|g);

    int n = (a|c|g|t);

    Int_matrix *bin2pos = new Int_matrix(n+1,"bin2pos");

    bin2pos->s(0,a);
    bin2pos->s(1,c);
    bin2pos->s(2,g);
    bin2pos->s(3,t);

    bin2pos->s(4,r);
    bin2pos->s(5,y);
    bin2pos->s(6,m);
    bin2pos->s(7,k);
    bin2pos->s(8,w);
    bin2pos->s(9,s);

    bin2pos->s(10,b);
    bin2pos->s(11,d);
    bin2pos->s(12,h);
    bin2pos->s(13,v);

    bin2pos->s(14,n);

    Int_matrix *pos2bin = new Int_matrix(n+1,"pos2bin");

    pos2bin->s(a,0);
    pos2bin->s(c,1);
    pos2bin->s(g,2);
    pos2bin->s(t,3);

    pos2bin->s(r,4);
    pos2bin->s(y,5);
    pos2bin->s(m,6);
    pos2bin->s(k,7);
    pos2bin->s(w,8);
    pos2bin->s(s,9);

    pos2bin->s(b,10);
    pos2bin->s(d,11);
    pos2bin->s(h,12);
    pos2bin->s(v,13);

    pos2bin->s(n,14);

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");
    mostcommon_table = new Int_matrix(char_fas,char_fas,"mostcommon_char"); // not used for DNA

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            int v = (pos2bin->g(i)&pos2bin->g(j));
            if(v>0)
            {
                parsimony_table->s(bin2pos->g(v),i,j);
                mostcommon_table->s(bin2pos->g(v),i,j);
            }
            else
            {
                parsimony_table->s(bin2pos->g((pos2bin->g(i)|pos2bin->g(j))),i,j);
                mostcommon_table->s(bin2pos->g((pos2bin->g(i)|pos2bin->g(j))),i,j);
            }
        }
    }


    // for situation where the parent's state has been updated and child's state may need to be changed.
    // if child's state is included in parent's state, use parsimony_table to get the minimum overlap;
    // if child's state is not included, change has happened between child and parent and child is not updated.

    child_parsimony_table = new Int_matrix(char_fas,char_fas,"child_parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            if( (pos2bin->g(i)&pos2bin->g(j)) >0)
            {
                int v = parsimony_table->g(i,j);
                child_parsimony_table->s(v,i,j);
            }
            else
            {
                child_parsimony_table->s(j,i,j);
            }
        }
    }



    delete bin2pos;
    delete pos2bin;

    if(Settings::noise>5)
    {
        stringstream ss;
        ss<<"\nModel_factory::define_dna_alphabet(). DNA parsimony table.\n\n  ";
        for(int i=0;i<15;i++)
            ss<<full_char_alphabet.at(i)<<" ";
        ss<<endl;

        for(int i=0;i<15;i++)
        {
            ss<<full_char_alphabet.at(i)<<" ";
            for(int j=0;j<15;j++)
            {
                ss<<full_char_alphabet.at(parsimony_table->g(i,j))<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        ss<<"\nModel_factory::define_dna_alphabet(). DNA child parsimony table.\n\n  ";
        for(int i=0;i<15;i++)
            ss<<full_char_alphabet.at(i)<<" ";
        ss<<endl;

        for(int i=0;i<15;i++)
        {
            ss<<full_char_alphabet.at(i)<<" ";
            for(int j=0;j<15;j++)
            {
                ss<<full_char_alphabet.at(child_parsimony_table->g(i,j))<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        Log_output::write_out(ss.str(),6);
    }


    char_ambiguity = new Db_matrix(4,char_fas);
    char_ambiguity->initialise(0);

}

/*******************************************/


void Model_factory::define_protein_alphabet()
{

    char_alphabet = this->protein_char_alphabet;
    full_char_alphabet = char_alphabet;

    char_as = char_alphabet.length();


    int count = 0;
    for(int i=0;i<char_as;i++)
    {
        Char_symbol letter;

        letter.index = i;
        letter.symbol = full_char_alphabet.at(i);
        letter.n_units = 1;
        letter.residues.push_back(full_char_alphabet.at(i));
        letter.first_residue = i;
        letter.second_residue = -1;
        char_symbols.push_back(letter);

        count++;
    }

    full_char_alphabet += "X";

    Char_symbol letter;

    letter.index = count;
    letter.symbol = 'X';
    letter.n_units = 20;
    for(int i=0;i<char_as;i++)
        letter.residues.push_back(full_char_alphabet.at(i));
    letter.first_residue = count;
    letter.second_residue = -1;
    char_symbols.push_back(letter);

    count++;

    for(int i=0;i<char_as-1;i++)
    {
        for(int j=i+1;j<char_as;j++)
        {
            Char_symbol letter;

            letter.index = count;
            letter.symbol = 'x';
            letter.n_units = 2;
            letter.residues.push_back(full_char_alphabet.at(i));
            letter.residues.push_back(full_char_alphabet.at(j));
            letter.first_residue = i;
            letter.second_residue = j;
            char_symbols.push_back(letter);

            full_char_alphabet += "X";

            count++;
        }
    }

    char_fas = full_char_alphabet.length();


//    this->print_char_alphabet();



    /************************************************************/

    double tmp_pi[20] = {0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466, 0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956};

    double tmp_q[400] = {-1.0644447077525, 0.024253680012, 0.0199296524112, 0.0421562148098, 0.019829882912, 0.0333710782038, 0.091898529865, 0.117944490096, 0.0077435982602, 0.00937017411, 0.034303854235, 0.056214349179, 0.0174255844392, 0.0080896843586, 0.065832507505, 0.234330242141, 0.129414648097, 0.0016275200247, 0.008491734537, 0.142217282556,
                         0.0477814374309, -0.9283507809291, 0.0248352939324, 0.0084029714104, 0.0101982061898, 0.11148814755, 0.0254969723473, 0.048674413647, 0.052213352795, 0.009062124214, 0.042903719239, 0.331941090612, 0.0133235035374, 0.0039473788809, 0.0310955230559, 0.085103118001, 0.0338262340451, 0.016744036728, 0.0134582713486, 0.0178549859644,
                         0.0441670615592, 0.027937434312, -1.38391347672512, 0.309721806842, 0.0051215097968, 0.056694964284, 0.0549932739622, 0.093704896008, 0.096657307877, 0.026861601976, 0.011338897352, 0.186830763486, 0.0038658446967, 0.00369569221099, 0.0089275113111, 0.276280123717, 0.123859441762, 0.00103458645453, 0.0383077812, 0.0139129779176,
                         0.0640178448442, 0.006477251488, 0.212232770148, -0.94297329251968, 0.00058492787022, 0.0226532677023, 0.358464938024, 0.0720614260512, 0.0227376245588, 0.001911353642, 0.0073109283823, 0.029764733853, 0.0020234831358, 0.00179593805976, 0.0194028221904, 0.074506504504, 0.0228715867982, 0.0018668150853, 0.0114891949562, 0.010799881226,
                         0.088970318416, 0.023225614652, 0.0103686978864, 0.00172817559999, -0.46436390434772, 0.00362939371299, 0.0012396736328, 0.0255311625132, 0.0060827096236, 0.00824576291, 0.033128997983, 0.00459221916954, 0.0076154533014, 0.015296664838, 0.0050066661924, 0.097857567114, 0.0312985388968, 0.010315697313, 0.0191832740086, 0.071047316584,
                         0.0787099366842, 0.133477006, 0.060339961416, 0.0351844479133, 0.00190795624962, -1.31468853172694, 0.317551411783, 0.0274774230936, 0.104910689643, 0.005521101322, 0.074957777201, 0.24159519414, 0.030136742202, 0.00384014619352, 0.0427139961732, 0.071524881773, 0.0523445036856, 0.0031035709083, 0.008032288082, 0.0213594972636,
                         0.137118971515, 0.019310611604, 0.0370254015012, 0.352205574616, 0.0004122601456, 0.200883241107, -1.17853838814971, 0.0472634621406, 0.0139264517825, 0.00617432607, 0.013298858967, 0.160308574698, 0.0061457688348, 0.00311812993141, 0.0312266801005, 0.0490058789081, 0.0501991141155, 0.0022522133463, 0.0069244312826, 0.0417384374836,
                         0.122727478488, 0.02570888938, 0.043997465064, 0.0493773258384, 0.0059212002572, 0.0121221828612, 0.0329610245313, -0.4741383934945, 0.006093410533, 0.0014757945466, 0.0052849306733, 0.0231712797588, 0.00339542007, 0.0019189431989, 0.011146518267, 0.093280508578, 0.0137786810791, 0.0048478037397, 0.0036545482168, 0.0132749884132,
                         0.0274570594166, 0.0939747598, 0.154649002326, 0.0530905054876, 0.0048071015816, 0.157714501491, 0.0330950244725, 0.020763831438, -0.9455247927498, 0.00669751654, 0.043058119558, 0.0552322503552, 0.0078818406807, 0.0261095183349, 0.0318601786938, 0.0514549945251, 0.0288777379989, 0.0037772913771, 0.136632497248, 0.0083910614248,
                         0.0167482050465, 0.008221840588, 0.0216647526984, 0.0022496876087, 0.003284932553, 0.0041839549677, 0.0073964135655, 0.00253502563518, 0.003376161347, -1.17498318692916, 0.27336615273, 0.0200868455952, 0.083031965142, 0.040717445093, 0.00457305166728, 0.022206797976, 0.088966278632, 0.0030567591897, 0.014821160614, 0.55449575628,
                         0.0344705408285, 0.021883589212, 0.0051413506032, 0.00483769259197, 0.0074197365386, 0.0319346789409, 0.0089563400907, 0.00510364337166, 0.0122025059606, 0.15368423202, -0.69175870029473, 0.015975776073, 0.094666495854, 0.081290001923, 0.0190303105564, 0.0239655313281, 0.0199280900994, 0.0095710687431, 0.0140609310556, 0.127636184504,
                         0.0785078337935, 0.23531264024, 0.117737663694, 0.0273733764605, 0.00142943173442, 0.14305227669, 0.150049162927, 0.0310993759044, 0.0217544113216, 0.015694841712, 0.022203558995, -1.07152398375532, 0.0182209045452, 0.0034141362684, 0.0254852873376, 0.067232846627, 0.084623394646, 0.0019781331795, 0.0047007809888, 0.0216539266904,
                         0.0774016821384, 0.030039999464, 0.0077483399574, 0.0059186573054, 0.0075393483596, 0.056754463806, 0.0182957528036, 0.01449413838, 0.0098736900133, 0.20634205636, 0.41846021018, 0.0579518322936, -1.2597234186325, 0.045758173097, 0.0078405461599, 0.0343352383995, 0.092502574724, 0.0074188949454, 0.0151127724254, 0.14593504782,
                         0.0182346531826, 0.004516408092, 0.00375891879174, 0.00266574034104, 0.007684890556, 0.00366990113448, 0.00471054498671, 0.0041568456258, 0.0165979167123, 0.05134827302, 0.18234669053, 0.0055103727096, 0.023220499701, -0.67999937105987, 0.0073881779164, 0.0379519766649, 0.0104882661681, 0.022005248076, 0.227669563576, 0.0460744832752,
                         0.124618565545, 0.029878490308, 0.0076255992414, 0.0241862096784, 0.0021123505512, 0.0342809801532, 0.0396167807095, 0.020277640926, 0.0170090221974, 0.0048431492208, 0.035849495396, 0.0345434792256, 0.0033413780883, 0.0062045996636, -0.5770185229931, 0.112151837712, 0.0485285253768, 0.0020054663895, 0.0076208498132, 0.0223241027972,
                         0.292004459041, 0.05383008268, 0.155350266162, 0.061138656376, 0.027178817748, 0.037788440247, 0.0409279829071, 0.111708930276, 0.0180832908897, 0.01548197904, 0.029719604451, 0.059989719918, 0.0096324810435, 0.0209811655989, 0.073828693968, -1.326554610767, 0.267114820854, 0.0075345000378, 0.0277605484806, 0.0165001710484,
                         0.183747304969, 0.024378648436, 0.079353827364, 0.0213842684566, 0.0099045924752, 0.0315100653768, 0.0477688308585, 0.0188010037494, 0.0115635053091, 0.07067118256, 0.028157755998, 0.086032427628, 0.029568433524, 0.0066065589057, 0.0363992375304, 0.304350756558, -1.1004826896859, 0.0015948784176, 0.0102700127816, 0.098419398788,
                         0.0098004742107, 0.05117989024, 0.00281118065298, 0.0074025714917, 0.013845044146, 0.0079236101097, 0.0090895272073, 0.0280544413194, 0.0064149020097, 0.010298201078, 0.057355623581, 0.008529242643, 0.0100576594062, 0.058786971516, 0.0063796049555, 0.0364094439818, 0.0067641119728, -0.44467569893618, 0.087670143938, 0.0259030544764,
                         0.0208543675065, 0.016776769076, 0.0424510884, 0.0185802165661, 0.0105002187974, 0.008363355651, 0.0113971362467, 0.0086252194872, 0.094633174672, 0.02036395922, 0.034364459162, 0.0082661793504, 0.0083556782799, 0.248050243532, 0.0098869347026, 0.0547101006747, 0.0177637255796, 0.035754572001, -0.6920103710931, 0.022312972188,
                         0.173776433679, 0.011074304228, 0.0076711383924, 0.0086899653085, 0.019349118692, 0.0110654786961, 0.0341810742559, 0.0155886497946, 0.0028916398054, 0.3790671258, 0.15520551106, 0.0189456434124, 0.040145332815, 0.0249765843548, 0.0144102052697, 0.0161795265281, 0.084699660521, 0.0052561618971, 0.011101848966, -1.034275403476 };

    Db_matrix *cPi = new Db_matrix(char_as,"pi_char");
    Db_matrix *cQ  = new Db_matrix(char_as,char_as,"Q_char");

    for(int j=0;j<char_as;j++)
    {
        cPi->s(tmp_pi[j],j);
        for(int i=0;i<char_as;i++)
        {
            cQ->s(tmp_q[j*char_as+i],j,i);
        }
    }


    /************************************************************/


    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            if(i==j)
            {
                parsimony_table->s(i,i,j);
            }
            else
            {
                Char_symbol *letter1 = &char_symbols.at(i);
                Char_symbol *letter2 = &char_symbols.at(j);

                if(letter1->index == char_as) // 'X'
                {
                    parsimony_table->s(j,i,j);
                }

                else if(letter2->index == char_as) // 'X'
                {
                    parsimony_table->s(i,i,j);
                }

                else if(letter1->n_units == 1 && letter2->n_units == 1)
                {
                    for(int k=char_as;k<char_fas;k++)
                    {
                        Char_symbol *letterX = &char_symbols.at(k);

                        if( ( letterX->first_residue == letter1->first_residue && letterX->second_residue == letter2->first_residue )
                            || ( letterX->first_residue == letter2->first_residue && letterX->second_residue == letter1->first_residue ) )
                        {
                            parsimony_table->s(k,i,j);
                        }
                    }
                }

                else if(letter1->n_units == 1 && letter2->n_units == 2 &&
                        ( letter1->first_residue == letter2->first_residue ||
                            letter1->first_residue == letter2->second_residue ) )
                {
                    parsimony_table->s(letter1->first_residue,i,j);
                }

                else if(letter2->n_units == 1 && letter1->n_units == 2 &&
                        ( letter2->first_residue == letter1->first_residue ||
                            letter2->first_residue == letter1->second_residue ) )
                {
                    parsimony_table->s(letter2->first_residue,i,j);
                }

                else
                {
                    // search for max score in Q matrix

                    float maxQ = -1; char maxl1 = 0; char maxl2 = 0;

                    int m = letter1->first_residue;
                    int n = letter2->first_residue;

                    if(cQ->g(m,n)>maxQ)
                    {
                        maxQ = cQ->g(m,n);
                        maxl1 = letter1->first_residue;
                        maxl2 = letter2->first_residue;
                    }

                    if(letter2->n_units == 2)
                    {

                        m = letter1->first_residue;
                        n = letter2->second_residue;

                        if(cQ->g(m,n)>maxQ)
                        {
                            maxQ = cQ->g(m,n);
                            maxl1 = letter1->first_residue;
                            maxl2 = letter2->second_residue;
                        }
                    }

                    if(letter1->n_units == 2)
                    {
                        m = letter1->second_residue;
                        n = letter2->first_residue;

                        if(cQ->g(m,n)>maxQ)
                        {
                            maxQ = cQ->g(m,n);
                            maxl1 = letter1->second_residue;
                            maxl2 = letter2->first_residue;
                        }
                    }

                    if(letter1->n_units == 2 && letter2->n_units == 2)
                    {
                        m = letter1->second_residue;
                        n = letter2->second_residue;

                        if(cQ->g(m,n)>maxQ)
                        {
                            maxQ = cQ->g(m,n);
                            maxl1 = letter1->second_residue;
                            maxl2 = letter2->second_residue;
                        }
                    }


                    // find the corresponding ambiguity character

                    for(int k=20;k<char_fas;k++)
                    {
                        Char_symbol *letterX = &char_symbols.at(k);

                        if( ( letterX->first_residue == maxl1 && letterX->second_residue == maxl2 )
                                || ( letterX->second_residue == maxl1 && letterX->first_residue == maxl2 ) )
                        {
                            parsimony_table->s(k,i,j);
                        }
                    }
                }

            }
        } // for(j)
    } // for(i)

    delete cPi;
    delete cQ;

    /*********************/

    // for situation where the parent's state has been updated and child's state may need to be changed.
    // if child's state is included in parent's state, use parsimony_table to get the minimum overlap;
    // if child's state is not included, change has happened between child and parent and child is not updated.

    child_parsimony_table = new Int_matrix(char_fas,char_fas,"child_parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            Char_symbol *parent = &char_symbols.at(i);
            Char_symbol *child = &char_symbols.at(j);

            // parent and child are identical
            if(i == j)
            {
                child_parsimony_table->s(j,i,j);
            }

            // one of them is X -> the other must have more or equal amount of information
            else if(parent->index == char_as) // 'X'
            {
                child_parsimony_table->s(j,i,j);
            }

            else if(child->index == char_as) // 'X'
            {
                child_parsimony_table->s(i,i,j);
            }

            // child is known: no change
            else if(child->n_units == 1)
            {
                child_parsimony_table->s(j,i,j);
            }

            // parent is known: only two comparisons
            else if(parent->n_units == 1 && parent->symbol != 'X')
            {
                if(parent->first_residue == child->first_residue || parent->first_residue == child->second_residue)
                    child_parsimony_table->s(i,i,j);
                else
                    child_parsimony_table->s(j,i,j);
            }

            // both are known: four comparisons
            else
            {
                int c = -1;
                if(parent->first_residue == child->first_residue || parent->first_residue == child->second_residue)
                    c = parent->first_residue;
                else if(parent->second_residue == child->first_residue || parent->second_residue == child->second_residue)
                    c = parent->second_residue;

                if(c>=0)
                {
                    for(int k=0;k<char_as;k++)
                    {
                        Char_symbol *letterX = &char_symbols.at(k);
                        if(letterX->first_residue == c || letterX->second_residue == c)
                        {
                            child_parsimony_table->s(k,i,j);
                        }
                    }
                }
                else
                {
                    child_parsimony_table->s(j,i,j);
//                    child_parsimony_table->s(-1,i,j);
                }

            }
        }
    }

    mostcommon_table = new Int_matrix(char_as,char_as,"mostcommon_char");
    for(int i=0;i<char_as;i++)
    {
        for(int j=0;j<char_as;j++)
        {
            mostcommon_table->s(j,i,j);
            if(tmp_pi[i]>tmp_pi[j])
                mostcommon_table->s(i,i,j);
        }
    }

//    this->print_int_matrix(mostcommon_table);

//    this->print_int_matrix(parsimony_table);

//    this->print_int_matrix(child_parsimony_table);


    if(Settings::noise>5)
    {
        stringstream ss;

        ss<<"\nModel_factory::define_protein_alphabet(). Protein parsimony table.\n\n  ";
        for(int i=0;i<char_fas;i++)
            ss<<full_char_alphabet.at(i)<<" ";
        ss<<endl;

        for(int i=0;i<char_fas;i++)
        {
            ss<<full_char_alphabet.at(i)<<" ";
            for(int j=0;j<char_fas;j++)
            {
                ss<<full_char_alphabet.at(parsimony_table->g(i,j))<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        Log_output::write_out(ss.str(),6);
    }

    char_ambiguity = new Db_matrix(20,char_fas);
    char_ambiguity->initialise(0);

}

/*
 * Definition of protein alpahbet including ambiguity characters.
 */
void Model_factory::define_protein_alphabet_groups()
{
    char_alphabet = "ARNDCQEGHILKMFPSTWYV";
    full_char_alphabet = "ARNDCQEGHILKMFPSTWYVZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZX";
    full_char_alphabet = "ARNDCQEGHILKMFPSTWYVabcdefghijklmnopqrstuvxyz12345X";

    char_as = char_alphabet.length();
    char_fas = full_char_alphabet.length();

    string ambiguity[] = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V",
                          "NG","HA","IV","ST","QE","ML","RK","ED","CV","FY","RKQ","AST","HML","NED","TIV",
                          "MLF","CIV","LFY","IML","NAST","HRKQ","ASTG","MLFY","LFYW","RKHSA","HRKQSTA",
                          "HRKQNEDSTA","HRKQNEDSTPA","HRKQNEDSTGPA","HRKQNEDSTGPACVIM","HRKQNEDSTGPACVIMLFYW"};

    for(int i=0;i<51;i++)
    {
        Char_symbol letter;

        letter.index = i;
        letter.symbol = full_char_alphabet.at(i);
        letter.n_units = ambiguity[i].length();
        for(int j=0;j<letter.n_units;j++)
            letter.residues.push_back(ambiguity[i].at(j));
        char_symbols.push_back(letter);
    }

    /*
     * Table for resolving the parental state from the parsimony ambiguity alphabet.
     */

    int table[] = {
        0, 44, 39, 46, 49, 45, 46, 41, 21, 49, 50, 44, 49, 50, 47, 31, 31, 50, 50, 49, 48,  0, 49, 31, 46, 50, 44, 46, 49, 50, 45,  0, 50, 46, 49, 50, 49, 50, 50,  0, 45,  0, 50, 50,  0,  0,  0,  0,  0,  0,  0,
       44,  1, 46, 46, 49, 30, 46, 48, 40, 49, 50, 26, 49, 50, 47, 44, 45, 50, 50, 49, 48, 44, 49, 45, 46, 50,  1, 46, 49, 50,  1, 45, 50, 46, 49, 50, 49, 50, 50, 46,  1, 48, 50, 50,  1,  1,  1,  1,  1,  1,  1,
       39, 46,  2, 33, 49, 46, 33, 20, 46, 49, 50, 46, 49, 50, 47, 39, 39, 50, 50, 49,  2, 46, 49, 39, 46, 50, 46, 33, 49, 50, 46, 39, 50,  2, 49, 50, 49, 50, 50,  2, 46, 48, 50, 50, 46, 46,  2,  2,  2,  2,  2,
       46, 46, 33,  3, 49, 46, 27, 48, 46, 49, 50, 46, 49, 50, 47, 46, 46, 50, 50, 49, 48, 46, 49, 46, 46, 50, 46,  3, 49, 50, 46, 46, 50,  3, 49, 50, 49, 50, 50, 46, 46, 48, 50, 50, 46, 46,  3,  3,  3,  3,  3,
       49, 49, 49, 49,  4, 49, 49, 49, 49, 36, 50, 49, 49, 50, 49, 49, 49, 50, 50, 28, 49, 49, 36, 49, 49, 50, 49, 49,  4, 50, 49, 49, 50, 49, 49, 50,  4, 50, 50, 49, 49, 49, 50, 50, 49, 49, 49, 49, 49,  4,  4,
       45, 30, 46, 46, 49,  5, 24, 48, 40, 49, 50, 30, 49, 50, 47, 45, 45, 50, 50, 49, 48, 45, 49, 45,  5, 50, 30, 46, 49, 50,  5, 45, 50, 46, 49, 50, 49, 50, 50, 46,  5, 48, 50, 50, 45,  5,  5,  5,  5,  5,  5,
       46, 46, 33, 27, 49, 24,  6, 48, 46, 49, 50, 46, 49, 50, 47, 46, 46, 50, 50, 49, 48, 46, 49, 46,  6, 50, 46,  6, 49, 50, 46, 46, 50,  6, 49, 50, 49, 50, 50, 46, 46, 48, 50, 50, 46, 46,  6,  6,  6,  6,  6,
       41, 48, 20, 48, 49, 48, 48,  7, 48, 49, 50, 48, 49, 50, 48, 41, 41, 50, 50, 49,  7, 48, 49, 41, 48, 50, 48, 48, 49, 50, 48, 41, 50, 48, 49, 50, 49, 50, 50, 48, 48,  7, 50, 50, 48, 48, 48, 48,  7,  7,  7,
       21, 40, 46, 46, 49, 40, 46, 48,  8, 49, 32, 40, 32, 50, 47, 44, 45, 50, 50, 49, 48,  8, 49, 45, 46, 32, 40, 46, 49, 50, 40, 45,  8, 46, 49, 50, 49, 50, 50, 46,  8, 48, 50, 50,  8,  8,  8,  8,  8,  8,  8,
       49, 49, 49, 49, 36, 49, 49, 49, 49,  9, 38, 49, 38, 50, 49, 49, 34, 50, 50, 22, 49, 49,  9, 49, 49, 38, 49, 49, 36, 50, 49, 49, 50, 49,  9, 50,  9, 50,  9, 49, 49, 49, 50, 50, 49, 49, 49, 49, 49,  9,  9,
       50, 50, 50, 50, 50, 50, 50, 50, 32, 38, 10, 50, 25, 35, 50, 50, 50, 43, 37, 50, 50, 50, 50, 50, 50, 10, 50, 50, 50, 37, 50, 50, 10, 50, 50, 10, 50, 10, 10, 50, 50, 50, 10, 10, 50, 50, 50, 50, 50, 50, 10,
       44, 26, 46, 46, 49, 30, 46, 48, 40, 49, 50, 11, 49, 50, 47, 44, 45, 50, 50, 49, 48, 44, 49, 45, 46, 50, 11, 46, 49, 50, 11, 45, 50, 46, 49, 50, 49, 50, 50, 46, 11, 48, 50, 50, 11, 11, 11, 11, 11, 11, 11,
       49, 49, 49, 49, 49, 49, 49, 49, 32, 38, 25, 49, 12, 35, 49, 49, 49, 50, 42, 49, 49, 49, 49, 49, 49, 12, 49, 49, 49, 42, 49, 49, 12, 49, 49, 12, 49, 42, 12, 49, 49, 49, 12, 50, 49, 49, 49, 49, 49, 12, 12,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 35, 50, 35, 13, 50, 50, 50, 43, 29, 50, 50, 50, 50, 50, 50, 35, 50, 50, 50, 13, 50, 50, 50, 50, 50, 13, 50, 13, 50, 50, 50, 50, 13, 13, 50, 50, 50, 50, 50, 50, 13,
       47, 47, 47, 47, 49, 47, 47, 48, 47, 49, 50, 47, 49, 50, 14, 47, 47, 50, 50, 49, 48, 47, 49, 47, 47, 50, 47, 47, 49, 50, 47, 47, 50, 47, 49, 50, 49, 50, 50, 47, 47, 48, 50, 50, 47, 47, 47, 14, 14, 14, 14,
       31, 44, 39, 46, 49, 45, 46, 41, 44, 49, 50, 44, 49, 50, 47, 15, 23, 50, 50, 49, 48, 44, 49, 15, 46, 50, 44, 46, 49, 50, 45, 15, 50, 46, 49, 50, 49, 50, 50, 15, 45, 15, 50, 50, 15, 15, 15, 15, 15, 15, 15,
       31, 45, 39, 46, 49, 45, 46, 41, 45, 34, 50, 45, 49, 50, 47, 23, 16, 50, 50, 34, 48, 45, 34, 16, 46, 50, 45, 46, 49, 50, 45, 16, 50, 46, 16, 50, 49, 50, 50, 16, 45, 16, 50, 50, 45, 16, 16, 16, 16, 16, 16,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 43, 50, 50, 43, 50, 50, 50, 17, 43, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 43, 50, 50, 50, 50, 50, 50, 50, 43, 50, 50, 50, 50, 50, 17, 50, 50, 50, 50, 50, 50, 17,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 37, 50, 42, 29, 50, 50, 50, 43, 18, 50, 50, 50, 50, 50, 50, 42, 50, 50, 50, 18, 50, 50, 50, 50, 50, 42, 50, 18, 50, 50, 50, 50, 18, 18, 50, 50, 50, 50, 50, 50, 18,
       49, 49, 49, 49, 28, 49, 49, 49, 49, 22, 50, 49, 49, 50, 49, 49, 34, 50, 50, 19, 49, 49, 19, 49, 49, 50, 49, 49, 19, 50, 49, 49, 50, 49, 19, 50, 19, 50, 50, 49, 49, 49, 50, 50, 49, 49, 49, 49, 49, 19, 19,
       48, 48,  2, 48, 49, 48, 48,  7, 48, 49, 50, 48, 49, 50, 48, 48, 48, 50, 50, 49, 20, 48, 49, 48, 48, 50, 48, 48, 49, 50, 48, 48, 50, 48, 49, 50, 49, 50, 50, 48, 48, 48, 50, 50, 48, 48, 48, 48, 20, 20, 20,
        0, 44, 46, 46, 49, 45, 46, 48,  8, 49, 50, 44, 49, 50, 47, 44, 45, 50, 50, 49, 48, 21, 49, 45, 46, 50, 44, 46, 49, 50, 45, 45, 50, 46, 49, 50, 49, 50, 50, 46, 45, 48, 50, 50, 21, 21, 21, 21, 21, 21, 21,
       49, 49, 49, 49, 36, 49, 49, 49, 49,  9, 50, 49, 49, 50, 49, 49, 34, 50, 50, 19, 49, 49, 22, 49, 49, 50, 49, 49, 36, 50, 49, 49, 50, 49, 22, 50, 22, 50, 50, 49, 49, 49, 50, 50, 49, 49, 49, 49, 49, 22, 22,
       31, 45, 39, 46, 49, 45, 46, 41, 45, 49, 50, 45, 49, 50, 47, 15, 16, 50, 50, 49, 48, 45, 49, 23, 46, 50, 45, 46, 49, 50, 45, 23, 50, 46, 49, 50, 49, 50, 50, 23, 45, 23, 50, 50, 45, 23, 23, 23, 23, 23, 23,
       46, 46, 46, 46, 49,  5,  6, 48, 46, 49, 50, 46, 49, 50, 47, 46, 46, 50, 50, 49, 48, 46, 49, 46, 24, 50, 46, 46, 49, 50, 46, 46, 50, 46, 49, 50, 49, 50, 50, 46, 46, 48, 50, 50, 46, 46, 24, 24, 24, 24, 24,
       50, 50, 50, 50, 50, 50, 50, 50, 32, 38, 10, 50, 12, 35, 50, 50, 50, 50, 42, 50, 50, 50, 50, 50, 50, 25, 50, 50, 50, 42, 50, 50, 25, 50, 50, 25, 50, 42, 25, 50, 50, 50, 25, 50, 50, 50, 50, 50, 50, 50, 25,
       44,  1, 46, 46, 49, 30, 46, 48, 40, 49, 50, 11, 49, 50, 47, 44, 45, 50, 50, 49, 48, 44, 49, 45, 46, 50, 26, 46, 49, 50, 26, 45, 50, 46, 49, 50, 49, 50, 50, 46, 26, 48, 50, 50, 26, 26, 26, 26, 26, 26, 26,
       46, 46, 33,  3, 49, 46,  6, 48, 46, 49, 50, 46, 49, 50, 47, 46, 46, 50, 50, 49, 48, 46, 49, 46, 46, 50, 46, 27, 49, 50, 46, 46, 50, 27, 49, 50, 49, 50, 50, 46, 46, 48, 50, 50, 46, 46, 27, 27, 27, 27, 27,
       49, 49, 49, 49,  4, 49, 49, 49, 49, 36, 50, 49, 49, 50, 49, 49, 49, 50, 50, 19, 49, 49, 36, 49, 49, 50, 49, 49, 28, 50, 49, 49, 50, 49, 49, 50, 28, 50, 50, 49, 49, 49, 50, 50, 49, 49, 49, 49, 49, 28, 28,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 37, 50, 42, 13, 50, 50, 50, 43, 18, 50, 50, 50, 50, 50, 50, 42, 50, 50, 50, 29, 50, 50, 50, 50, 50, 42, 50, 29, 50, 50, 50, 50, 29, 29, 50, 50, 50, 50, 50, 50, 29,
       45,  1, 46, 46, 49,  5, 46, 48, 40, 49, 50, 11, 49, 50, 47, 45, 45, 50, 50, 49, 48, 45, 49, 45, 46, 50, 26, 46, 49, 50, 30, 45, 50, 46, 49, 50, 49, 50, 50, 46, 30, 48, 50, 50, 45, 30, 30, 30, 30, 30, 30,
        0, 45, 39, 46, 49, 45, 46, 41, 45, 49, 50, 45, 49, 50, 47, 15, 16, 50, 50, 49, 48, 45, 49, 23, 46, 50, 45, 46, 49, 50, 45, 31, 50, 46, 49, 50, 49, 50, 50, 31, 45, 31, 50, 50, 45, 31, 31, 31, 31, 31, 31,
       50, 50, 50, 50, 50, 50, 50, 50,  8, 50, 10, 50, 12, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 25, 50, 50, 50, 50, 50, 50, 32, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 32,
       46, 46,  2,  3, 49, 46,  6, 48, 46, 49, 50, 46, 49, 50, 47, 46, 46, 50, 50, 49, 48, 46, 49, 46, 46, 50, 46, 27, 49, 50, 46, 46, 50, 33, 49, 50, 49, 50, 50, 46, 46, 48, 50, 50, 46, 46, 33, 33, 33, 33, 33,
       49, 49, 49, 49, 49, 49, 49, 49, 49,  9, 50, 49, 49, 50, 49, 49, 16, 50, 50, 19, 49, 49, 22, 49, 49, 50, 49, 49, 49, 50, 49, 49, 50, 49, 34, 50, 49, 50, 50, 49, 49, 49, 50, 50, 49, 49, 49, 49, 49, 34, 34,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 10, 50, 12, 13, 50, 50, 50, 50, 42, 50, 50, 50, 50, 50, 50, 25, 50, 50, 50, 42, 50, 50, 50, 50, 50, 35, 50, 42, 50, 50, 50, 50, 35, 50, 50, 50, 50, 50, 50, 50, 35,
       49, 49, 49, 49,  4, 49, 49, 49, 49,  9, 50, 49, 49, 50, 49, 49, 49, 50, 50, 19, 49, 49, 22, 49, 49, 50, 49, 49, 28, 50, 49, 49, 50, 49, 49, 50, 36, 50, 50, 49, 49, 49, 50, 50, 49, 49, 49, 49, 49, 36, 36,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 10, 50, 42, 13, 50, 50, 50, 43, 18, 50, 50, 50, 50, 50, 50, 42, 50, 50, 50, 29, 50, 50, 50, 50, 50, 42, 50, 37, 50, 50, 50, 50, 37, 37, 50, 50, 50, 50, 50, 50, 37,
       50, 50, 50, 50, 50, 50, 50, 50, 50,  9, 10, 50, 12, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 25, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 38, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 38,
        0, 46,  2, 46, 49, 46, 46, 48, 46, 49, 50, 46, 49, 50, 47, 15, 16, 50, 50, 49, 48, 46, 49, 23, 46, 50, 46, 46, 49, 50, 46, 31, 50, 46, 49, 50, 49, 50, 50, 39, 46, 48, 50, 50, 46, 46, 39, 39, 39, 39, 39,
       45,  1, 46, 46, 49,  5, 46, 48,  8, 49, 50, 11, 49, 50, 47, 45, 45, 50, 50, 49, 48, 45, 49, 45, 46, 50, 26, 46, 49, 50, 30, 45, 50, 46, 49, 50, 49, 50, 50, 46, 40, 48, 50, 50, 45, 40, 40, 40, 40, 40, 40,
        0, 48, 48, 48, 49, 48, 48,  7, 48, 49, 50, 48, 49, 50, 48, 15, 16, 50, 50, 49, 48, 48, 49, 23, 48, 50, 48, 48, 49, 50, 48, 31, 50, 48, 49, 50, 49, 50, 50, 48, 48, 41, 50, 50, 48, 48, 48, 48, 41, 41, 41,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 10, 50, 12, 13, 50, 50, 50, 50, 18, 50, 50, 50, 50, 50, 50, 25, 50, 50, 50, 29, 50, 50, 50, 50, 50, 35, 50, 37, 50, 50, 50, 50, 42, 50, 50, 50, 50, 50, 50, 50, 42,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 10, 50, 50, 13, 50, 50, 50, 17, 18, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 29, 50, 50, 50, 50, 50, 50, 50, 37, 50, 50, 50, 50, 50, 43, 50, 50, 50, 50, 50, 50, 43,
        0,  1, 46, 46, 49, 45, 46, 48,  8, 49, 50, 11, 49, 50, 47, 15, 45, 50, 50, 49, 48, 21, 49, 45, 46, 50, 26, 46, 49, 50, 45, 45, 50, 46, 49, 50, 49, 50, 50, 46, 45, 48, 50, 50, 44, 44, 44, 44, 44, 44, 44,
        0,  1, 46, 46, 49,  5, 46, 48,  8, 49, 50, 11, 49, 50, 47, 15, 16, 50, 50, 49, 48, 21, 49, 23, 46, 50, 26, 46, 49, 50, 30, 31, 50, 46, 49, 50, 49, 50, 50, 46, 40, 48, 50, 50, 44, 45, 45, 45, 45, 45, 45,
        0,  1,  2,  3, 49,  5,  6, 48,  8, 49, 50, 11, 49, 50, 47, 15, 16, 50, 50, 49, 48, 21, 49, 23, 24, 50, 26, 27, 49, 50, 30, 31, 50, 33, 49, 50, 49, 50, 50, 39, 40, 48, 50, 50, 44, 45, 46, 46, 46, 46, 46,
        0,  1,  2,  3, 49,  5,  6, 48,  8, 49, 50, 11, 49, 50, 14, 15, 16, 50, 50, 49, 48, 21, 49, 23, 24, 50, 26, 27, 49, 50, 30, 31, 50, 33, 49, 50, 49, 50, 50, 39, 40, 48, 50, 50, 44, 45, 46, 47, 47, 47, 47,
        0,  1,  2,  3, 49,  5,  6,  7,  8, 49, 50, 11, 49, 50, 14, 15, 16, 50, 50, 49, 20, 21, 49, 23, 24, 50, 26, 27, 49, 50, 30, 31, 50, 33, 49, 50, 49, 50, 50, 39, 40, 41, 50, 50, 44, 45, 46, 47, 48, 48, 48,
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 50, 11, 12, 50, 14, 15, 16, 50, 50, 19, 20, 21, 22, 23, 24, 50, 26, 27, 28, 50, 30, 31, 50, 33, 34, 50, 36, 50, 50, 39, 40, 41, 50, 50, 44, 45, 46, 47, 48, 49, 49,
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    };

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");
    mostcommon_table = new Int_matrix(char_fas,char_fas,"mostcommon_char"); // not used for groups

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            parsimony_table->s(table[i*char_fas+j],i,j);
            mostcommon_table->s(table[i*char_fas+j],i,j);
        }
    }

    // for situation where the parent's state has been updated and child's state may need to be changed.
    // if child's state is included in parent's state, use parsimony_table to get the minimum overlap;
    // if child's state is not included, change has happened between child and parent and child is not updated.

    child_parsimony_table = new Int_matrix(char_fas,char_fas,"child_parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            string parent_amb = ambiguity[i];
            string child_amb = ambiguity[j];

            bool all_included = true;


            for(int k=0;k<(int)parent_amb.size();k++)
            {
                if(child_amb.find(parent_amb.at(k))==string::npos)
                {
                    all_included = false;
                    break;
                }
            }


            if( all_included )
            {
                int v = parsimony_table->g(i,j);
                child_parsimony_table->s(v,i,j);
            }
            else
            {
                child_parsimony_table->s(j,i,j);
            }
        }
    }

    mostcommon_table = new Int_matrix(char_as,char_as,"mostcommon_char");

    if(Settings::noise>5)
    {
        stringstream ss;

        ss<<"\nModel_factory::define_protein_alphabet(). Protein parsimony table.\n\n  ";
        for(int i=0;i<char_fas;i++)
            ss<<full_char_alphabet.at(i)<<" ";
        ss<<endl;

        for(int i=0;i<char_fas;i++)
        {
            ss<<full_char_alphabet.at(i)<<" ";
            for(int j=0;j<char_fas;j++)
            {
                ss<<full_char_alphabet.at(parsimony_table->g(i,j))<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        Log_output::write_out(ss.str(),6);
    }

    char_ambiguity = new Db_matrix(20,char_fas);
    char_ambiguity->initialise(0);

}

/*******************************************/

/*
 * Definition of protein alpahbet including ambiguity characters.
 */
void Model_factory::define_codon_alphabet()
{
    string full_alpha = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTNNN";

    char_alphabet = "";
    full_char_alphabet = "";

    char_as  = 61;

    int count = 0;

    for(int i=0;i<61;i++)
    {
        Codon_symbol codon;

        codon.symbol = full_alpha.substr(i*3,3);
        codon.index = i;
        codon.n_units = 1;
        codon.codons.push_back(codon.symbol);
        codon.first_codon = i;
        codon.second_codon = -1;
        codon_symbols.push_back(codon);

        count++;
    }


    Codon_symbol codon;

    codon.symbol = "NNN";
    codon.index = 61;
    codon.n_units = 1;
    for(int i=0;i<61;i++)
        codon.codons.push_back(full_alpha.substr(i*3,3));
    codon.first_codon = 61;
    codon.second_codon = -1;
    codon_symbols.push_back(codon);

    count++;

    for(int i=0;i<60;i++)
    {
        for(int j=i+1;j<61;j++)
        {
            Codon_symbol codon;

            codon.index = count;
            codon.symbol = "nnn";
            codon.n_units = 2;
            codon.codons.push_back(full_alpha.substr(i*3,3));
            codon.codons.push_back(full_alpha.substr(j*3,3));
            codon.first_codon = i;
            codon.second_codon = j;
            codon_symbols.push_back(codon);

            count++;
        }
    }

    char_fas = count;

    double tmp_pi[61] = {0.024709, 0.026990, 0.037138, 0.018760, 0.013945, 0.024015, 0.009522, 0.014230, 0.010544, 0.015952, 0.008214, 0.008979, 0.007064, 0.029436, 0.020724, 0.019431, 0.015416, 0.013015, 0.026398, 0.008551, 0.011465, 0.012756, 0.006648, 0.009507, 0.005177, 0.012128, 0.007718, 0.007807, 0.005328, 0.017807, 0.031660, 0.011119, 0.028694, 0.027083, 0.035815, 0.023890, 0.013035, 0.027437, 0.008784, 0.019503, 0.015756, 0.021140, 0.009781, 0.015576, 0.007598, 0.019661, 0.026815, 0.015815, 0.022021, 0.013288, 0.009289, 0.015669, 0.007124, 0.012150, 0.019672, 0.015708, 0.014291, 0.007146, 0.026531, 0.012766, 0.015811};

    double tmp_q[3721] = {-0.996852, 0.019663, 0.434175, 0.025324, 0.017695, 0.012761, 0.005660, 0.011280, 0.061068, 0.009776, 0.021005, 0.009536, 0.005768, 0.002711, 0.011898, 0.003740, 0.042628, 0.005906, 0.026732, 0.007654, 0.007352, 0.002504, 0.002391, 0.003500, 0.012882, 0.015172, 0.009349, 0.014707, 0.001452, 0.001182, 0.001249, 0.003244, 0.061214, 0.009708, 0.018164, 0.012120, 0.009128, 0.004818, 0.003400, 0.009158, 0.004921, 0.003252, 0.000831, 0.002467, 0.003219, 0.004301, 0.002816, 0.003567, 0.001069, 0.002186, 0.008633, 0.004091, 0.002042, 0.006228, 0.000140, 0.001428, 0.000110, 0.002805, 0.000981, 0.002652, 0.001441,
                          0.018001, -0.833627, 0.030232, 0.299934, 0.007529, 0.029937, 0.005116, 0.010141, 0.003332, 0.059604, 0.002731, 0.016622, 0.001258, 0.003052, 0.006465, 0.001080, 0.008318, 0.033977, 0.017907, 0.012141, 0.002650, 0.004096, 0.001733, 0.002844, 0.002169, 0.005276, 0.001991, 0.003428, 0.000548, 0.001018, 0.001301, 0.000545, 0.009888, 0.070669, 0.013486, 0.023938, 0.002984, 0.013093, 0.003870, 0.004338, 0.008120, 0.023487, 0.007724, 0.006311, 0.001142, 0.002033, 0.001953, 0.001510, 0.007353, 0.001702, 0.005085, 0.012603, 0.004073, 0.005464, 0.001329, 0.000791, 0.000847, 0.000644, 0.002019, 0.000763, 0.001431,
                          0.288865, 0.021971, -0.843432, 0.014913, 0.007281, 0.018314, 0.011510, 0.007363, 0.023323, 0.011365, 0.053669, 0.004764, 0.001533, 0.002715, 0.014706, 0.001548, 0.016566, 0.007204, 0.058416, 0.005632, 0.003355, 0.003986, 0.003516, 0.002606, 0.007632, 0.031929, 0.024684, 0.011513, 0.001164, 0.002390, 0.004846, 0.001680, 0.018442, 0.013350, 0.042174, 0.006639, 0.004700, 0.008692, 0.007474, 0.007676, 0.004622, 0.006910, 0.002998, 0.001889, 0.002176, 0.004544, 0.005439, 0.002552, 0.002017, 0.002332, 0.003315, 0.006478, 0.006181, 0.004989, 0.000223, 0.000560, 0.000073, 0.000814, 0.000961, 0.003282, 0.000970,
                          0.033355, 0.431526, 0.029523, -0.987099, 0.014903, 0.010645, 0.007384, 0.023703, 0.006355, 0.020501, 0.004741, 0.036958, 0.001072, 0.001118, 0.006785, 0.003785, 0.020131, 0.010786, 0.014616, 0.024978, 0.004115, 0.002128, 0.001443, 0.003259, 0.002156, 0.001891, 0.002442, 0.003918, 0.000564, 0.000919, 0.001150, 0.002896, 0.020626, 0.031842, 0.011045, 0.072138, 0.004404, 0.004479, 0.003079, 0.007530, 0.012646, 0.007896, 0.006704, 0.014683, 0.001177, 0.001220, 0.002259, 0.004108, 0.003652, 0.006112, 0.009177, 0.006368, 0.004857, 0.011702, 0.000900, 0.000509, 0.001259, 0.000831, 0.001006, 0.002572, 0.002573,
                          0.031354, 0.014573, 0.019391, 0.020049, -1.243500, 0.259288, 0.147445, 0.260513, 0.015627, 0.024450, 0.004176, 0.020780, 0.021256, 0.005638, 0.021888, 0.008172, 0.016821, 0.003028, 0.012152, 0.006987, 0.007839, 0.003930, 0.004451, 0.003892, 0.002448, 0.002018, 0.003008, 0.002370, 0.002550, 0.002950, 0.004344, 0.003024, 0.017739, 0.005225, 0.010136, 0.006291, 0.021145, 0.009424, 0.008924, 0.012785, 0.003689, 0.002751, 0.001829, 0.002262, 0.011715, 0.010648, 0.016584, 0.011067, 0.002333, 0.003536, 0.048425, 0.021486, 0.017250, 0.020995, 0.000984, 0.000529, 0.001355, 0.003769, 0.004284, 0.006642, 0.003291,
                          0.013129, 0.033646, 0.028321, 0.008316, 0.150556, -1.097510, 0.139816, 0.197266, 0.002217, 0.059113, 0.002412, 0.014527, 0.003885, 0.034135, 0.020238, 0.004692, 0.007230, 0.009949, 0.019241, 0.003572, 0.001964, 0.007404, 0.002603, 0.002921, 0.001748, 0.009255, 0.003360, 0.001426, 0.000805, 0.009008, 0.006028, 0.001764, 0.009818, 0.013657, 0.014062, 0.004799, 0.006814, 0.044137, 0.007396, 0.009147, 0.001874, 0.005363, 0.002093, 0.001265, 0.005544, 0.032109, 0.015770, 0.007693, 0.006001, 0.001725, 0.008444, 0.051019, 0.016328, 0.014919, 0.002173, 0.001355, 0.000560, 0.001103, 0.007393, 0.002331, 0.002065,
                          0.014687, 0.014502, 0.044894, 0.014547, 0.215929, 0.352630, -1.406200, 0.189339, 0.004280, 0.030649, 0.007156, 0.016698, 0.011924, 0.013213, 0.028368, 0.008267, 0.009388, 0.005630, 0.030042, 0.005121, 0.003732, 0.005286, 0.006377, 0.003060, 0.002635, 0.004339, 0.005676, 0.002939, 0.001577, 0.006201, 0.010942, 0.002943, 0.010931, 0.012754, 0.021891, 0.008278, 0.009608, 0.026715, 0.021540, 0.010194, 0.002835, 0.004400, 0.003624, 0.001051, 0.006111, 0.017671, 0.037522, 0.009466, 0.006702, 0.003924, 0.014147, 0.026283, 0.036125, 0.017230, 0.003900, 0.002086, 0.001308, 0.001115, 0.007715, 0.005218, 0.002888,
                          0.019586, 0.019234, 0.019218, 0.031249, 0.255292, 0.332923, 0.126698, -1.266280, 0.002834, 0.026528, 0.003170, 0.034494, 0.004782, 0.010430, 0.010621, 0.022796, 0.012602, 0.005524, 0.009239, 0.007390, 0.002554, 0.002426, 0.001468, 0.006049, 0.001136, 0.002621, 0.001335, 0.002962, 0.001345, 0.002380, 0.002930, 0.002783, 0.016675, 0.005733, 0.010980, 0.012087, 0.009933, 0.016212, 0.005285, 0.027461, 0.003338, 0.002825, 0.002368, 0.003573, 0.007188, 0.008487, 0.010058, 0.023212, 0.002478, 0.002207, 0.023887, 0.020905, 0.010822, 0.041072, 0.000922, 0.001716, 0.002281, 0.001789, 0.002582, 0.003275, 0.002328,
                          0.143109, 0.008529, 0.082151, 0.011307, 0.020667, 0.005050, 0.003865, 0.003825, -1.172030, 0.013413, 0.230080, 0.008486, 0.004219, 0.000918, 0.008596, 0.001573, 0.019807, 0.004337, 0.017706, 0.006355, 0.004002, 0.001226, 0.001383, 0.002244, 0.119639, 0.086404, 0.077589, 0.187931, 0.002005, 0.001361, 0.002177, 0.002336, 0.011718, 0.002232, 0.005356, 0.004616, 0.004753, 0.003025, 0.001352, 0.004934, 0.009323, 0.002279, 0.003855, 0.001965, 0.003069, 0.001752, 0.002632, 0.002507, 0.001718, 0.002002, 0.005187, 0.001598, 0.001824, 0.002826, 0.000740, 0.001169, 0.000452, 0.002036, 0.000829, 0.003273, 0.000720,
                          0.015143, 0.100850, 0.026460, 0.024110, 0.021374, 0.088994, 0.018295, 0.023664, 0.008865, -1.134840, 0.011222, 0.184047, 0.001266, 0.005081, 0.008771, 0.002019, 0.010616, 0.015081, 0.025476, 0.005128, 0.005533, 0.009222, 0.004831, 0.005342, 0.004001, 0.018797, 0.006563, 0.001963, 0.001404, 0.003401, 0.005569, 0.000883, 0.013753, 0.044440, 0.024604, 0.010764, 0.006736, 0.035682, 0.008081, 0.007624, 0.010821, 0.047683, 0.014155, 0.010143, 0.001711, 0.006961, 0.005633, 0.002857, 0.007246, 0.002259, 0.023880, 0.083006, 0.042204, 0.035529, 0.012829, 0.001691, 0.002213, 0.000597, 0.003560, 0.001680, 0.002524,
                          0.063186, 0.008973, 0.242659, 0.010827, 0.007089, 0.007053, 0.008295, 0.005492, 0.295346, 0.021795, -1.428430, 0.013076, 0.002855, 0.002348, 0.026156, 0.001463, 0.013660, 0.008049, 0.052049, 0.005814, 0.002194, 0.005230, 0.002925, 0.002948, 0.064463, 0.136976, 0.196339, 0.077079, 0.000977, 0.002183, 0.008981, 0.002248, 0.004789, 0.005236, 0.016869, 0.002400, 0.004462, 0.006344, 0.005525, 0.004035, 0.004823, 0.004205, 0.015277, 0.002281, 0.001683, 0.001982, 0.010400, 0.001254, 0.003289, 0.003215, 0.002211, 0.004172, 0.003565, 0.002322, 0.001597, 0.009974, 0.000339, 0.000816, 0.002412, 0.003565, 0.000659,
                          0.026242, 0.049966, 0.019705, 0.077218, 0.032273, 0.038856, 0.017708, 0.054666, 0.009965, 0.326984, 0.011962, -1.277050, 0.002464, 0.002193, 0.006083, 0.007151, 0.017106, 0.008477, 0.026180, 0.010608, 0.004151, 0.002994, 0.001776, 0.005592, 0.003719, 0.006082, 0.002849, 0.004285, 0.001048, 0.001740, 0.003428, 0.002061, 0.022995, 0.022963, 0.018742, 0.030636, 0.017795, 0.016301, 0.008316, 0.021800, 0.015141, 0.022784, 0.014023, 0.019830, 0.001705, 0.002613, 0.003840, 0.002359, 0.002096, 0.006557, 0.051613, 0.050652, 0.044803, 0.075464, 0.003223, 0.001183, 0.006296, 0.000676, 0.002201, 0.001503, 0.003405,
                          0.020177, 0.004806, 0.008062, 0.002846, 0.041962, 0.013209, 0.016074, 0.009634, 0.006297, 0.002858, 0.003319, 0.003132, -1.366290, 0.307296, 0.071081, 0.317317, 0.004254, 0.002661, 0.004990, 0.002242, 0.001804, 0.001521, 0.002119, 0.002147, 0.002004, 0.001584, 0.001030, 0.001793, 0.012765, 0.020674, 0.042996, 0.013209, 0.010283, 0.001041, 0.004226, 0.000884, 0.008055, 0.005836, 0.003125, 0.005231, 0.000923, 0.001197, 0.000587, 0.000291, 0.079194, 0.042651, 0.087968, 0.059847, 0.003415, 0.004866, 0.003778, 0.002968, 0.001696, 0.002722, 0.001694, 0.001694, 0.000923, 0.037165, 0.011791, 0.026013, 0.010364,
                          0.002275, 0.002798, 0.003425, 0.000712, 0.002671, 0.027849, 0.004274, 0.005042, 0.000329, 0.002754, 0.000655, 0.000669, 0.073742, -0.887525, 0.036537, 0.277085, 0.000588, 0.001591, 0.002204, 0.000938, 0.000462, 0.000916, 0.000432, 0.000442, 0.000310, 0.001433, 0.000481, 0.000205, 0.005975, 0.051463, 0.049826, 0.011136, 0.000537, 0.001757, 0.001098, 0.000491, 0.001488, 0.008330, 0.001383, 0.001952, 0.000316, 0.000749, 0.000322, 0.000261, 0.011818, 0.116456, 0.077458, 0.032484, 0.003766, 0.001830, 0.000435, 0.001592, 0.000552, 0.000764, 0.003157, 0.001140, 0.001163, 0.005466, 0.023741, 0.011331, 0.006466,
                          0.014186, 0.008420, 0.026354, 0.006142, 0.014728, 0.023453, 0.013034, 0.007293, 0.004374, 0.006752, 0.010367, 0.002635, 0.024228, 0.051896, -0.819425, 0.034272, 0.009840, 0.004428, 0.033136, 0.005408, 0.001793, 0.002788, 0.002853, 0.001316, 0.001677, 0.003878, 0.003750, 0.001634, 0.015012, 0.042873, 0.100407, 0.022120, 0.005320, 0.001841, 0.010676, 0.001606, 0.007612, 0.017861, 0.009238, 0.012917, 0.003152, 0.004325, 0.003224, 0.000943, 0.009435, 0.016533, 0.059619, 0.012872, 0.009247, 0.003839, 0.004118, 0.006914, 0.003795, 0.002646, 0.005194, 0.007425, 0.001392, 0.017283, 0.022359, 0.043380, 0.013614,
                          0.004756, 0.001500, 0.002959, 0.003654, 0.005865, 0.005799, 0.004051, 0.016694, 0.000854, 0.001658, 0.000618, 0.003304, 0.115351, 0.419742, 0.036551, -1.067930, 0.002964, 0.001118, 0.001734, 0.002224, 0.000596, 0.000262, 0.000268, 0.001307, 0.000385, 0.000430, 0.000447, 0.000391, 0.009769, 0.015524, 0.031405, 0.039917, 0.002721, 0.000739, 0.001379, 0.001514, 0.004182, 0.003910, 0.001435, 0.006983, 0.000798, 0.000572, 0.000571, 0.000505, 0.021872, 0.045457, 0.058521, 0.112609, 0.001667, 0.002998, 0.000872, 0.000815, 0.000619, 0.002007, 0.001237, 0.000746, 0.001347, 0.014550, 0.007522, 0.019678, 0.018014,
                          0.068322, 0.014562, 0.039907, 0.024496, 0.015215, 0.011262, 0.005798, 0.011632, 0.013547, 0.010984, 0.007278, 0.009963, 0.001949, 0.001123, 0.013227, 0.003736, -1.042850, 0.017118, 0.362704, 0.020271, 0.012857, 0.002382, 0.004217, 0.004555, 0.018466, 0.010309, 0.008920, 0.008461, 0.003882, 0.002895, 0.006793, 0.007226, 0.096006, 0.010655, 0.041474, 0.019065, 0.012615, 0.006118, 0.003430, 0.011218, 0.007635, 0.002822, 0.002564, 0.003833, 0.003724, 0.003495, 0.006193, 0.004666, 0.003032, 0.005180, 0.010133, 0.010580, 0.004960, 0.012746, 0.000311, 0.004037, 0.000749, 0.005155, 0.001154, 0.003769, 0.001467,
                          0.011212, 0.070462, 0.020558, 0.015546, 0.003244, 0.018358, 0.004119, 0.006040, 0.003514, 0.018484, 0.005080, 0.005848, 0.001444, 0.003598, 0.007051, 0.001670, 0.020276, -0.982175, 0.086164, 0.296699, 0.002026, 0.007498, 0.005144, 0.002620, 0.005918, 0.041083, 0.014976, 0.006844, 0.001322, 0.012436, 0.006408, 0.002790, 0.009192, 0.032040, 0.018704, 0.010491, 0.002675, 0.007635, 0.002781, 0.002644, 0.002334, 0.007338, 0.003462, 0.001320, 0.001756, 0.007376, 0.005956, 0.001372, 0.069267, 0.022290, 0.003399, 0.013588, 0.005789, 0.006074, 0.001508, 0.005598, 0.001004, 0.000925, 0.020558, 0.001278, 0.005392,
                          0.025021, 0.018309, 0.082182, 0.010387, 0.006419, 0.017504, 0.010836, 0.004980, 0.007072, 0.015395, 0.016195, 0.008904, 0.001335, 0.002457, 0.026013, 0.001277, 0.211815, 0.042481, -0.912145, 0.019229, 0.005688, 0.005581, 0.005953, 0.005359, 0.007026, 0.023328, 0.028657, 0.005375, 0.002017, 0.007659, 0.022625, 0.003488, 0.025534, 0.012897, 0.104029, 0.008146, 0.004521, 0.013093, 0.009729, 0.005362, 0.003817, 0.006040, 0.006087, 0.000923, 0.001709, 0.002517, 0.008074, 0.001461, 0.005150, 0.002767, 0.004550, 0.009790, 0.008290, 0.004397, 0.000303, 0.003554, 0.000442, 0.001441, 0.002415, 0.003586, 0.000954,
                          0.022117, 0.038319, 0.024460, 0.054797, 0.011395, 0.010031, 0.005702, 0.012297, 0.007836, 0.009566, 0.005585, 0.011138, 0.001852, 0.003228, 0.013105, 0.005053, 0.036545, 0.451575, 0.059362, -1.154830, 0.005724, 0.007304, 0.002793, 0.011668, 0.007607, 0.013909, 0.012681, 0.014874, 0.001300, 0.006561, 0.004202, 0.006791, 0.012162, 0.014662, 0.014606, 0.017232, 0.004931, 0.006277, 0.006123, 0.004896, 0.008044, 0.005020, 0.005979, 0.002388, 0.003028, 0.003561, 0.008095, 0.002696, 0.030126, 0.057219, 0.008048, 0.007923, 0.003978, 0.013302, 0.001946, 0.006667, 0.000645, 0.003330, 0.007293, 0.002892, 0.012386,
                          0.015844, 0.006239, 0.010867, 0.006734, 0.009535, 0.004114, 0.003100, 0.003170, 0.003680, 0.007699, 0.001572, 0.003251, 0.001111, 0.001186, 0.003241, 0.001010, 0.017288, 0.002300, 0.013097, 0.004270, -0.959559, 0.223948, 0.145555, 0.273040, 0.003185, 0.001680, 0.003125, 0.003955, 0.002936, 0.003838, 0.004289, 0.001923, 0.018380, 0.008443, 0.014106, 0.011581, 0.018027, 0.005523, 0.006333, 0.009839, 0.003453, 0.003361, 0.002370, 0.001994, 0.001710, 0.001503, 0.002242, 0.003604, 0.001639, 0.001006, 0.026870, 0.006940, 0.005643, 0.011351, 0.000348, 0.001150, 0.000191, 0.001656, 0.000677, 0.001652, 0.001186,
                          0.004851, 0.008666, 0.011606, 0.003130, 0.004296, 0.013940, 0.003946, 0.002706, 0.001014, 0.011533, 0.003368, 0.002108, 0.000842, 0.002114, 0.004529, 0.000399, 0.002879, 0.007651, 0.011550, 0.004897, 0.201287, -0.935790, 0.169605, 0.242298, 0.001797, 0.008043, 0.004252, 0.002579, 0.001284, 0.014390, 0.006850, 0.002471, 0.005843, 0.011779, 0.014689, 0.003705, 0.005502, 0.032614, 0.005974, 0.005071, 0.002405, 0.007984, 0.002819, 0.001219, 0.001962, 0.006785, 0.004305, 0.002786, 0.002103, 0.001358, 0.005064, 0.029821, 0.006204, 0.006575, 0.000827, 0.000599, 0.000124, 0.000385, 0.004301, 0.001462, 0.000644,
                          0.008887, 0.007035, 0.019642, 0.004072, 0.009336, 0.009404, 0.009134, 0.003141, 0.002194, 0.011591, 0.003614, 0.002398, 0.002252, 0.001915, 0.008893, 0.000783, 0.009780, 0.010071, 0.023640, 0.003593, 0.251028, 0.325436, -1.182370, 0.186111, 0.004204, 0.006686, 0.005666, 0.003038, 0.002577, 0.004869, 0.011430, 0.002942, 0.014556, 0.014253, 0.023689, 0.006993, 0.013193, 0.017819, 0.025836, 0.008778, 0.003664, 0.010686, 0.006029, 0.001110, 0.002483, 0.005086, 0.008847, 0.002354, 0.001142, 0.001789, 0.009471, 0.012101, 0.015888, 0.008013, 0.001294, 0.002566, 0.000513, 0.001184, 0.001734, 0.003184, 0.002748,
                          0.009096, 0.008075, 0.010179, 0.006431, 0.005709, 0.007379, 0.003065, 0.009054, 0.002489, 0.008963, 0.002547, 0.005281, 0.001595, 0.001369, 0.002868, 0.002672, 0.007386, 0.003586, 0.014879, 0.010495, 0.329258, 0.325081, 0.130133, -1.102160, 0.002969, 0.002771, 0.001888, 0.005911, 0.002293, 0.002767, 0.003912, 0.009884, 0.009055, 0.009660, 0.007452, 0.012225, 0.009565, 0.008649, 0.003323, 0.021290, 0.003185, 0.004371, 0.003693, 0.002225, 0.000986, 0.003033, 0.003503, 0.003594, 0.002146, 0.003193, 0.007029, 0.006857, 0.007335, 0.022330, 0.000201, 0.001501, 0.000337, 0.001061, 0.001215, 0.001115, 0.002049,
                          0.061482, 0.011309, 0.054753, 0.007813, 0.006592, 0.008107, 0.004847, 0.003122, 0.243662, 0.012327, 0.102276, 0.006450, 0.002734, 0.001762, 0.006711, 0.001444, 0.054989, 0.014878, 0.035828, 0.012564, 0.007053, 0.004428, 0.005399, 0.005453, -1.458750, 0.251947, 0.223508, 0.179170, 0.004224, 0.002538, 0.006210, 0.003035, 0.006631, 0.005195, 0.011259, 0.003910, 0.005948, 0.006714, 0.002751, 0.003063, 0.005696, 0.005783, 0.003366, 0.003023, 0.003245, 0.002399, 0.003680, 0.001829, 0.002609, 0.001321, 0.007931, 0.004985, 0.004540, 0.005820, 0.000648, 0.001916, 0.000737, 0.001275, 0.001801, 0.003068, 0.000989,
                          0.030911, 0.011743, 0.097776, 0.002925, 0.002320, 0.018328, 0.003406, 0.003076, 0.075119, 0.024724, 0.092770, 0.004503, 0.000922, 0.003479, 0.006626, 0.000689, 0.013104, 0.044088, 0.050777, 0.009807, 0.001589, 0.008459, 0.003665, 0.002172, 0.107549, -1.141450, 0.182320, 0.204540, 0.001315, 0.010476, 0.007677, 0.002957, 0.004510, 0.005214, 0.012346, 0.002033, 0.001404, 0.010192, 0.002825, 0.002443, 0.001313, 0.005931, 0.001950, 0.001221, 0.001125, 0.006989, 0.003346, 0.001499, 0.010181, 0.001788, 0.001883, 0.011372, 0.004751, 0.002943, 0.002711, 0.004452, 0.000215, 0.000722, 0.003107, 0.002172, 0.001000,
                          0.029932, 0.006962, 0.118786, 0.005935, 0.005434, 0.010457, 0.007003, 0.002461, 0.106003, 0.013565, 0.208965, 0.003315, 0.000942, 0.001835, 0.010070, 0.001125, 0.017819, 0.025256, 0.098023, 0.014051, 0.004642, 0.007029, 0.004881, 0.002326, 0.149932, 0.286510, -1.403230, 0.121475, 0.001491, 0.003138, 0.015265, 0.002298, 0.003852, 0.003249, 0.020813, 0.001939, 0.002406, 0.007198, 0.003845, 0.004052, 0.002357, 0.005332, 0.006612, 0.001242, 0.001155, 0.001502, 0.009330, 0.001263, 0.002561, 0.001966, 0.002016, 0.003848, 0.005643, 0.001572, 0.003441, 0.009626, 0.001129, 0.001244, 0.002761, 0.003215, 0.001132,
                          0.046550, 0.011851, 0.054770, 0.009416, 0.004233, 0.004387, 0.003585, 0.005400, 0.253828, 0.004012, 0.081101, 0.004928, 0.001622, 0.000774, 0.004338, 0.000973, 0.016709, 0.011410, 0.018176, 0.016293, 0.005808, 0.004214, 0.002587, 0.007199, 0.118820, 0.317765, 0.120091, -1.214070, 0.001018, 0.001176, 0.003451, 0.004755, 0.005100, 0.002322, 0.006418, 0.003142, 0.002528, 0.003315, 0.001071, 0.003160, 0.001590, 0.002880, 0.001780, 0.001881, 0.001542, 0.002053, 0.002506, 0.001905, 0.004527, 0.003029, 0.002849, 0.002085, 0.002827, 0.005805, 0.000532, 0.000421, 0.001925, 0.001248, 0.001258, 0.002375, 0.000762,
                          0.006733, 0.002776, 0.008117, 0.001988, 0.006675, 0.003629, 0.002818, 0.003591, 0.003968, 0.004203, 0.001506, 0.001766, 0.016925, 0.033010, 0.058397, 0.035632, 0.011234, 0.003230, 0.009993, 0.002086, 0.006319, 0.003074, 0.003216, 0.004092, 0.004104, 0.002995, 0.002160, 0.001492, -1.584530, 0.199432, 0.444164, 0.184670, 0.002556, 0.002069, 0.003553, 0.001793, 0.005501, 0.005072, 0.003964, 0.004594, 0.001764, 0.004516, 0.000774, 0.000726, 0.016601, 0.014883, 0.029983, 0.012714, 0.005580, 0.005128, 0.005481, 0.003655, 0.001979, 0.002564, 0.002159, 0.001855, 0.002083, 0.180994, 0.020807, 0.150360, 0.016826,
                          0.001640, 0.001543, 0.004984, 0.000968, 0.002310, 0.012149, 0.003316, 0.001902, 0.000806, 0.003047, 0.001007, 0.000877, 0.008201, 0.085071, 0.049895, 0.016940, 0.002507, 0.009090, 0.011354, 0.003151, 0.002471, 0.010308, 0.001818, 0.001478, 0.000738, 0.007135, 0.001360, 0.000516, 0.059667, -1.170780, 0.335384, 0.150181, 0.001178, 0.001656, 0.003341, 0.000602, 0.002488, 0.018559, 0.003134, 0.003252, 0.000579, 0.002579, 0.000989, 0.000651, 0.004389, 0.045858, 0.027239, 0.006564, 0.006871, 0.001753, 0.000914, 0.005262, 0.001523, 0.001822, 0.004500, 0.002840, 0.001328, 0.044211, 0.059985, 0.113983, 0.010919,
                          0.000975, 0.001109, 0.005684, 0.000681, 0.001913, 0.004572, 0.003291, 0.001317, 0.000725, 0.002806, 0.002330, 0.000972, 0.009593, 0.046325, 0.065722, 0.019274, 0.003308, 0.002634, 0.018865, 0.001135, 0.001553, 0.002760, 0.002400, 0.001175, 0.001015, 0.002941, 0.003721, 0.000851, 0.074741, 0.188631, -0.904017, 0.096500, 0.001086, 0.001079, 0.002661, 0.000256, 0.002255, 0.005230, 0.007045, 0.002148, 0.000437, 0.001095, 0.000855, 0.000231, 0.003592, 0.013058, 0.054897, 0.005046, 0.004292, 0.001680, 0.001265, 0.001862, 0.002565, 0.000824, 0.002028, 0.003398, 0.000498, 0.037640, 0.021434, 0.141955, 0.014083,
                          0.007209, 0.001323, 0.005610, 0.004886, 0.003793, 0.003810, 0.002521, 0.003561, 0.002215, 0.001267, 0.001660, 0.001665, 0.008391, 0.029481, 0.041228, 0.069758, 0.010019, 0.003265, 0.008282, 0.005223, 0.001983, 0.002835, 0.001759, 0.008452, 0.001413, 0.003226, 0.001595, 0.003338, 0.088483, 0.240513, 0.274777, -1.339180, 0.003745, 0.000851, 0.002341, 0.001435, 0.005900, 0.005917, 0.002263, 0.009201, 0.000952, 0.000804, 0.000579, 0.000868, 0.006795, 0.012082, 0.015077, 0.026331, 0.003221, 0.006827, 0.001106, 0.002110, 0.000827, 0.002982, 0.000486, 0.001516, 0.001824, 0.149304, 0.016458, 0.179996, 0.033836,
                          0.052712, 0.009301, 0.023870, 0.013485, 0.008621, 0.008217, 0.003627, 0.008269, 0.004306, 0.007646, 0.001371, 0.007195, 0.002531, 0.000551, 0.003842, 0.001842, 0.051581, 0.004169, 0.023491, 0.003625, 0.007344, 0.002598, 0.003372, 0.003000, 0.001196, 0.001906, 0.001036, 0.001388, 0.000475, 0.000731, 0.001198, 0.001451, -0.909708, 0.051581, 0.414485, 0.075046, 0.012436, 0.008018, 0.005118, 0.013648, 0.011244, 0.003260, 0.002304, 0.002853, 0.004717, 0.001663, 0.004077, 0.003948, 0.001111, 0.002447, 0.008217, 0.003565, 0.003744, 0.005991, 0.000044, 0.001035, 0.000062, 0.000913, 0.000668, 0.001035, 0.000532,
                          0.008857, 0.070427, 0.018306, 0.022056, 0.002690, 0.012110, 0.004484, 0.003012, 0.000869, 0.026175, 0.001588, 0.007613, 0.000272, 0.001910, 0.001409, 0.000530, 0.006065, 0.015397, 0.012571, 0.004629, 0.003574, 0.005548, 0.003498, 0.003391, 0.000993, 0.002335, 0.000926, 0.000669, 0.000407, 0.001088, 0.001261, 0.000349, 0.054649, -0.887159, 0.113244, 0.366986, 0.003726, 0.013919, 0.003567, 0.004327, 0.006999, 0.025744, 0.007086, 0.005164, 0.000556, 0.003034, 0.001970, 0.001121, 0.002863, 0.001857, 0.002944, 0.009274, 0.003594, 0.005240, 0.000088, 0.001469, 0.000123, 0.000460, 0.001223, 0.000518, 0.000405,
                          0.012531, 0.010163, 0.043732, 0.005785, 0.003946, 0.009429, 0.005820, 0.004362, 0.001577, 0.010958, 0.003869, 0.004699, 0.000833, 0.000902, 0.006177, 0.000748, 0.017852, 0.006797, 0.076676, 0.003487, 0.004516, 0.005231, 0.004397, 0.001978, 0.001627, 0.004181, 0.004485, 0.001399, 0.000529, 0.001661, 0.002352, 0.000727, 0.332074, 0.085634, -0.827623, 0.041741, 0.005780, 0.017439, 0.013747, 0.007609, 0.002938, 0.004955, 0.006328, 0.001732, 0.002173, 0.003341, 0.010932, 0.002150, 0.001869, 0.000733, 0.003117, 0.006317, 0.005270, 0.004059, 0.000168, 0.000740, 0.000058, 0.000648, 0.000494, 0.001538, 0.000613,
                          0.012536, 0.027045, 0.010321, 0.056647, 0.003672, 0.004824, 0.003299, 0.007199, 0.002037, 0.007188, 0.000825, 0.011514, 0.000261, 0.000605, 0.001393, 0.001232, 0.012303, 0.005715, 0.009001, 0.006168, 0.005558, 0.001978, 0.001946, 0.004865, 0.000847, 0.001032, 0.000627, 0.001027, 0.000400, 0.000449, 0.000340, 0.000668, 0.090138, 0.416040, 0.062578, -0.870644, 0.003111, 0.005763, 0.002636, 0.008120, 0.011860, 0.007154, 0.004118, 0.010285, 0.000834, 0.000876, 0.001336, 0.004188, 0.001246, 0.002364, 0.007991, 0.007123, 0.003352, 0.011668, 0.000187, 0.000394, 0.000066, 0.000267, 0.000650, 0.001134, 0.001641,
                          0.017303, 0.006179, 0.013392, 0.006338, 0.022620, 0.012553, 0.007019, 0.010843, 0.003845, 0.008243, 0.002812, 0.012257, 0.004365, 0.003361, 0.012101, 0.006234, 0.014920, 0.002671, 0.009156, 0.003235, 0.015856, 0.005384, 0.006728, 0.006976, 0.002362, 0.001306, 0.001424, 0.001514, 0.002248, 0.003399, 0.005476, 0.005033, 0.027375, 0.007742, 0.015882, 0.005701, -1.281480, 0.271180, 0.121716, 0.349548, 0.025036, 0.012072, 0.009808, 0.009767, 0.019693, 0.013435, 0.033049, 0.016723, 0.001508, 0.001516, 0.033088, 0.018297, 0.011147, 0.019128, 0.002102, 0.002326, 0.004203, 0.003712, 0.002389, 0.004174, 0.002013,
                          0.004339, 0.012881, 0.011766, 0.003063, 0.004790, 0.038633, 0.009271, 0.008408, 0.001163, 0.020746, 0.001899, 0.005334, 0.001503, 0.008937, 0.013491, 0.002769, 0.003438, 0.003622, 0.012598, 0.001956, 0.002308, 0.015163, 0.004317, 0.002997, 0.001267, 0.004505, 0.002025, 0.000943, 0.000985, 0.012045, 0.006035, 0.002398, 0.008386, 0.013740, 0.022765, 0.005018, 0.128837, -1.030450, 0.130039, 0.243982, 0.006334, 0.030818, 0.005660, 0.007301, 0.003441, 0.044243, 0.025613, 0.007859, 0.005323, 0.001728, 0.009658, 0.050009, 0.013015, 0.016277, 0.008300, 0.001580, 0.001481, 0.000944, 0.008794, 0.001789, 0.001926,
                          0.009566, 0.011891, 0.031600, 0.006577, 0.014167, 0.020223, 0.023350, 0.008562, 0.001623, 0.014675, 0.005167, 0.008501, 0.002513, 0.004636, 0.021795, 0.003175, 0.006020, 0.004121, 0.029241, 0.005961, 0.008266, 0.008675, 0.019554, 0.003597, 0.001622, 0.003900, 0.003378, 0.000952, 0.002404, 0.006353, 0.025396, 0.002864, 0.016719, 0.010999, 0.056056, 0.007170, 0.180631, 0.406194, -1.458140, 0.207610, 0.008643, 0.024944, 0.009536, 0.008580, 0.008847, 0.025380, 0.041698, 0.007592, 0.003135, 0.001237, 0.015863, 0.024816, 0.022260, 0.023773, 0.006593, 0.000907, 0.002668, 0.002343, 0.004267, 0.006127, 0.003222,
                          0.011602, 0.006004, 0.014616, 0.007243, 0.009141, 0.011263, 0.004977, 0.020036, 0.002668, 0.006236, 0.001699, 0.010036, 0.001894, 0.002947, 0.013725, 0.006957, 0.008868, 0.001765, 0.007258, 0.002147, 0.005784, 0.003316, 0.002992, 0.010379, 0.000813, 0.001519, 0.001603, 0.001265, 0.001255, 0.002969, 0.003487, 0.005246, 0.020080, 0.006009, 0.013972, 0.009947, 0.233624, 0.343229, 0.093501, -1.138310, 0.008293, 0.010152, 0.005295, 0.016253, 0.006961, 0.015340, 0.017163, 0.031158, 0.002042, 0.001846, 0.019846, 0.024005, 0.009984, 0.038159, 0.003614, 0.000576, 0.004391, 0.001681, 0.002354, 0.003321, 0.003799,
                          0.007718, 0.013910, 0.010894, 0.015057, 0.003265, 0.002856, 0.001713, 0.003015, 0.006239, 0.010956, 0.002514, 0.008628, 0.000414, 0.000590, 0.004146, 0.000985, 0.007471, 0.001928, 0.006395, 0.004366, 0.002513, 0.001947, 0.001546, 0.001922, 0.001872, 0.001011, 0.001155, 0.000788, 0.000597, 0.000654, 0.000879, 0.000672, 0.020477, 0.012030, 0.006679, 0.017983, 0.020712, 0.011030, 0.004818, 0.010265, -1.001990, 0.242481, 0.193001, 0.290581, 0.000570, 0.001342, 0.001968, 0.001121, 0.001111, 0.000979, 0.007586, 0.003922, 0.003928, 0.011724, 0.001627, 0.001603, 0.001521, 0.000366, 0.001458, 0.000658, 0.001838,
                          0.003801, 0.029987, 0.012139, 0.007007, 0.001814, 0.006092, 0.001982, 0.001902, 0.001136, 0.035981, 0.001634, 0.009677, 0.000400, 0.001044, 0.004240, 0.000525, 0.002058, 0.004518, 0.007542, 0.002031, 0.001823, 0.004818, 0.003360, 0.001966, 0.001416, 0.003403, 0.001947, 0.001064, 0.001138, 0.002173, 0.001640, 0.000423, 0.004426, 0.032982, 0.008394, 0.008085, 0.007444, 0.039997, 0.010364, 0.009366, 0.180725, -0.893499, 0.162152, 0.220989, 0.000952, 0.004233, 0.004363, 0.001619, 0.001751, 0.001882, 0.002320, 0.010545, 0.004139, 0.005933, 0.001842, 0.002790, 0.000292, 0.000318, 0.003491, 0.000432, 0.000992,
                          0.002099, 0.021315, 0.011383, 0.012858, 0.002607, 0.005139, 0.003528, 0.003445, 0.004156, 0.023086, 0.012829, 0.012873, 0.000424, 0.000968, 0.006831, 0.001134, 0.004042, 0.004606, 0.016430, 0.005227, 0.002779, 0.003677, 0.004098, 0.003590, 0.001782, 0.002418, 0.005217, 0.001421, 0.000421, 0.001800, 0.002767, 0.000658, 0.006759, 0.019620, 0.023171, 0.010059, 0.013071, 0.015877, 0.008564, 0.010558, 0.310904, 0.350466, -1.194700, 0.187362, 0.001220, 0.002572, 0.012307, 0.001980, 0.001040, 0.001711, 0.003135, 0.006574, 0.004036, 0.006286, 0.000900, 0.005809, 0.000218, 0.000354, 0.001905, 0.001128, 0.001507,
                          0.003913, 0.010936, 0.004503, 0.017684, 0.002025, 0.001950, 0.000642, 0.003265, 0.001330, 0.010388, 0.001203, 0.011430, 0.000132, 0.000493, 0.001255, 0.000630, 0.003793, 0.001103, 0.001564, 0.001311, 0.001467, 0.000998, 0.000474, 0.001358, 0.001005, 0.000951, 0.000616, 0.000943, 0.000248, 0.000744, 0.000470, 0.000620, 0.005256, 0.008979, 0.003983, 0.015775, 0.008173, 0.012861, 0.004838, 0.020351, 0.293935, 0.299925, 0.117652, -0.910018, 0.000359, 0.002136, 0.001832, 0.001882, 0.000719, 0.001709, 0.002387, 0.002737, 0.002236, 0.008580, 0.000484, 0.000704, 0.001019, 0.000314, 0.000642, 0.000455, 0.000651,
                          0.010467, 0.004056, 0.010637, 0.002905, 0.021502, 0.017523, 0.007659, 0.013462, 0.004259, 0.003591, 0.001820, 0.002015, 0.073626, 0.045786, 0.025733, 0.055936, 0.007556, 0.003008, 0.005937, 0.003408, 0.002580, 0.003293, 0.002172, 0.001234, 0.002211, 0.001796, 0.001173, 0.001584, 0.011640, 0.010286, 0.014968, 0.009944, 0.017814, 0.001984, 0.010243, 0.002624, 0.033787, 0.012424, 0.010227, 0.017869, 0.001183, 0.002650, 0.001570, 0.000735, -1.493460, 0.203898, 0.345114, 0.368061, 0.002358, 0.003131, 0.004999, 0.004287, 0.003022, 0.003004, 0.003421, 0.002632, 0.011550, 0.013059, 0.008267, 0.014380, 0.005397,
                          0.005405, 0.002791, 0.008583, 0.001164, 0.007552, 0.039220, 0.008558, 0.006143, 0.000940, 0.005648, 0.000828, 0.001193, 0.015323, 0.174351, 0.017427, 0.044926, 0.002741, 0.004883, 0.003379, 0.001549, 0.000876, 0.004402, 0.001720, 0.001467, 0.000632, 0.004311, 0.000589, 0.000815, 0.004033, 0.041533, 0.021028, 0.006833, 0.002427, 0.004179, 0.006087, 0.001065, 0.008907, 0.061739, 0.011338, 0.015217, 0.001076, 0.004552, 0.001280, 0.001692, 0.078794, -1.212990, 0.303670, 0.203144, 0.007099, 0.002857, 0.001471, 0.004969, 0.001776, 0.001700, 0.006424, 0.001449, 0.003064, 0.002814, 0.019319, 0.009426, 0.004614,
                          0.002595, 0.001966, 0.007534, 0.001581, 0.008624, 0.014124, 0.013324, 0.005337, 0.001035, 0.003351, 0.003186, 0.001286, 0.023173, 0.085029, 0.046076, 0.042407, 0.003561, 0.002891, 0.007949, 0.002581, 0.000959, 0.002048, 0.002193, 0.001242, 0.000710, 0.001513, 0.002685, 0.000730, 0.005957, 0.018089, 0.064818, 0.006252, 0.004363, 0.001990, 0.014601, 0.001190, 0.016066, 0.026207, 0.013659, 0.012483, 0.001156, 0.003439, 0.004489, 0.001064, 0.097787, 0.222658, -0.992739, 0.131449, 0.004976, 0.002213, 0.000945, 0.003502, 0.002068, 0.002101, 0.002046, 0.001658, 0.000851, 0.003054, 0.011043, 0.014682, 0.006192,
                          0.005572, 0.002577, 0.005992, 0.004873, 0.009758, 0.011682, 0.005699, 0.020885, 0.001672, 0.002882, 0.000651, 0.001340, 0.026730, 0.060459, 0.016867, 0.138357, 0.004548, 0.001129, 0.002439, 0.001458, 0.002612, 0.002247, 0.000990, 0.002160, 0.000599, 0.001149, 0.000616, 0.000940, 0.004283, 0.007391, 0.010101, 0.018512, 0.007163, 0.001919, 0.004868, 0.006327, 0.013783, 0.013634, 0.004217, 0.038423, 0.001116, 0.002164, 0.001224, 0.001854, 0.176821, 0.252543, 0.222870, -1.191070, 0.001942, 0.002855, 0.002745, 0.002012, 0.002085, 0.004058, 0.003119, 0.001003, 0.004932, 0.008134, 0.007219, 0.010875, 0.013991,
                          0.001200, 0.009012, 0.003402, 0.003111, 0.001477, 0.006545, 0.002898, 0.001601, 0.000823, 0.005249, 0.001227, 0.000855, 0.001095, 0.005034, 0.008702, 0.001471, 0.002123, 0.040938, 0.006174, 0.011699, 0.000853, 0.001218, 0.000345, 0.000926, 0.000613, 0.005607, 0.000897, 0.001605, 0.001350, 0.005556, 0.006171, 0.001627, 0.001448, 0.003521, 0.003039, 0.001352, 0.000893, 0.006632, 0.001250, 0.001809, 0.000795, 0.001680, 0.000462, 0.000509, 0.000814, 0.006338, 0.006059, 0.001394, -0.644620, 0.302981, 0.001357, 0.007620, 0.001670, 0.002468, 0.004368, 0.015225, 0.000389, 0.000956, 0.100380, 0.002562, 0.025244,
                          0.004065, 0.003457, 0.006518, 0.008629, 0.003710, 0.003118, 0.002812, 0.002363, 0.001588, 0.002712, 0.001987, 0.004431, 0.002587, 0.004053, 0.005988, 0.004384, 0.006010, 0.021832, 0.005496, 0.036822, 0.000868, 0.001304, 0.000895, 0.002285, 0.000515, 0.001632, 0.001142, 0.001780, 0.002056, 0.002350, 0.004004, 0.005713, 0.005285, 0.003786, 0.001976, 0.004250, 0.001487, 0.003567, 0.000818, 0.002709, 0.001161, 0.002995, 0.001259, 0.002003, 0.001790, 0.004227, 0.004466, 0.003398, 0.502104, -0.889242, 0.003720, 0.004493, 0.002479, 0.009767, 0.001267, 0.018544, 0.004176, 0.003871, 0.047826, 0.002770, 0.085944,
                          0.022964, 0.014775, 0.013255, 0.018534, 0.072694, 0.021831, 0.014502, 0.036591, 0.005888, 0.041008, 0.001955, 0.049888, 0.002873, 0.001379, 0.009188, 0.001824, 0.016817, 0.004762, 0.012932, 0.007409, 0.033163, 0.006954, 0.006778, 0.007195, 0.004420, 0.002458, 0.001675, 0.002394, 0.003144, 0.001753, 0.004313, 0.001324, 0.025381, 0.008583, 0.012017, 0.020550, 0.046431, 0.028525, 0.014999, 0.041667, 0.012867, 0.005280, 0.003301, 0.004002, 0.004089, 0.003114, 0.002728, 0.004673, 0.003216, 0.005322, -1.246680, 0.192909, 0.098011, 0.239222, 0.001093, 0.001864, 0.002915, 0.005421, 0.003352, 0.006357, 0.002154,
                          0.006451, 0.021710, 0.015353, 0.007624, 0.019121, 0.078194, 0.015972, 0.018985, 0.001075, 0.084504, 0.002187, 0.029025, 0.001338, 0.002990, 0.009145, 0.001011, 0.010410, 0.011286, 0.016494, 0.004324, 0.005078, 0.024276, 0.005134, 0.004161, 0.001647, 0.008802, 0.001895, 0.001039, 0.001243, 0.005980, 0.003763, 0.001498, 0.006529, 0.016030, 0.014439, 0.010861, 0.015221, 0.087565, 0.013911, 0.029879, 0.003944, 0.014226, 0.004104, 0.002721, 0.002079, 0.006236, 0.005994, 0.002031, 0.010709, 0.003810, 0.114363, -1.153740, 0.131575, 0.202385, 0.006684, 0.004617, 0.002342, 0.000649, 0.014746, 0.001673, 0.002707,
                          0.007082, 0.015432, 0.032219, 0.012789, 0.033765, 0.055041, 0.048284, 0.021615, 0.002700, 0.094500, 0.004110, 0.056467, 0.001681, 0.002281, 0.011038, 0.001689, 0.010734, 0.010575, 0.030719, 0.004774, 0.009082, 0.011108, 0.014826, 0.009788, 0.003299, 0.008087, 0.006113, 0.003097, 0.001480, 0.003806, 0.011401, 0.001291, 0.015080, 0.013661, 0.026493, 0.011241, 0.020395, 0.050122, 0.027444, 0.027333, 0.008686, 0.012283, 0.005541, 0.004889, 0.003222, 0.004900, 0.007782, 0.004629, 0.005161, 0.004624, 0.127796, 0.289390, -1.480360, 0.239634, 0.005493, 0.005746, 0.001772, 0.001551, 0.005466, 0.004601, 0.004546,
                          0.012665, 0.012137, 0.015248, 0.018067, 0.024096, 0.029488, 0.013503, 0.048101, 0.002453, 0.046646, 0.001570, 0.055765, 0.001583, 0.001851, 0.004513, 0.003210, 0.016172, 0.006506, 0.009553, 0.009362, 0.010711, 0.006902, 0.004384, 0.017472, 0.002480, 0.002938, 0.000999, 0.003730, 0.001124, 0.002670, 0.002147, 0.002729, 0.014149, 0.011680, 0.011963, 0.022941, 0.020521, 0.036754, 0.017185, 0.061251, 0.015203, 0.010322, 0.005060, 0.011000, 0.001879, 0.002751, 0.004636, 0.005281, 0.004473, 0.010682, 0.182890, 0.260996, 0.140506, -1.271040, 0.001663, 0.001754, 0.004489, 0.001936, 0.006448, 0.001867, 0.009983,
                          0.000176, 0.001824, 0.000421, 0.000858, 0.000698, 0.002653, 0.001888, 0.000667, 0.000397, 0.010403, 0.000667, 0.001471, 0.000608, 0.004724, 0.005471, 0.001222, 0.000244, 0.000998, 0.000407, 0.000846, 0.000203, 0.000536, 0.000437, 0.000097, 0.000171, 0.001672, 0.001350, 0.000211, 0.000585, 0.004073, 0.003265, 0.000274, 0.000064, 0.000121, 0.000306, 0.000227, 0.001393, 0.011576, 0.002944, 0.003583, 0.001303, 0.001980, 0.000448, 0.000383, 0.001321, 0.006420, 0.002789, 0.002507, 0.004890, 0.000856, 0.000516, 0.005324, 0.001989, 0.001027, -0.466903, 0.000873, 0.357340, 0.000823, 0.005093, 0.000286, 0.001006,
                          0.002247, 0.001360, 0.001324, 0.000608, 0.000469, 0.002072, 0.001264, 0.001555, 0.000784, 0.001717, 0.005216, 0.000676, 0.000762, 0.002136, 0.009796, 0.000923, 0.003962, 0.004638, 0.005973, 0.003629, 0.000839, 0.000486, 0.001086, 0.000908, 0.000632, 0.003437, 0.004729, 0.000209, 0.000629, 0.003220, 0.006850, 0.001073, 0.001891, 0.002532, 0.001688, 0.000599, 0.001931, 0.002759, 0.000507, 0.000715, 0.001608, 0.003755, 0.003617, 0.000698, 0.001273, 0.001814, 0.002831, 0.001010, 0.021345, 0.015687, 0.001102, 0.004606, 0.002606, 0.001357, 0.001093, -0.196689, 0.001299, 0.002416, 0.019373, 0.007356, 0.014010,
                          0.000190, 0.001600, 0.000190, 0.001653, 0.001322, 0.000941, 0.000872, 0.002271, 0.000333, 0.002470, 0.000195, 0.003956, 0.000456, 0.002396, 0.002018, 0.001832, 0.000808, 0.000914, 0.000817, 0.000386, 0.000153, 0.000111, 0.000239, 0.000224, 0.000267, 0.000182, 0.000610, 0.001051, 0.000777, 0.001654, 0.001104, 0.001419, 0.000124, 0.000234, 0.000144, 0.000110, 0.003833, 0.002843, 0.001640, 0.005992, 0.001677, 0.000433, 0.000149, 0.001110, 0.006141, 0.004216, 0.001596, 0.005458, 0.000599, 0.003883, 0.001895, 0.002568, 0.000883, 0.003817, 0.491875, 0.001428, -0.582349, 0.001114, 0.002527, 0.001361, 0.001287,
                          0.009698, 0.002434, 0.004231, 0.002181, 0.007355, 0.003706, 0.001486, 0.003563, 0.003005, 0.001332, 0.000938, 0.000850, 0.036737, 0.022514, 0.050124, 0.039565, 0.011123, 0.001684, 0.005324, 0.003985, 0.002657, 0.000688, 0.001101, 0.001412, 0.000924, 0.001225, 0.001344, 0.001364, 0.134940, 0.110170, 0.166768, 0.232319, 0.003665, 0.001742, 0.003248, 0.000893, 0.006771, 0.003626, 0.002881, 0.004588, 0.000807, 0.000941, 0.000485, 0.000684, 0.013885, 0.007741, 0.011459, 0.018003, 0.002946, 0.007198, 0.007047, 0.001422, 0.001547, 0.003293, 0.002265, 0.005311, 0.002229, -1.350170, 0.022013, 0.315940, 0.030797,
                          0.000913, 0.002054, 0.001345, 0.000711, 0.002252, 0.006692, 0.002769, 0.001385, 0.000329, 0.002140, 0.000747, 0.000745, 0.003139, 0.026339, 0.017464, 0.005509, 0.000670, 0.010085, 0.002403, 0.002351, 0.000292, 0.002068, 0.000435, 0.000435, 0.000351, 0.001420, 0.000803, 0.000370, 0.004178, 0.040259, 0.025578, 0.006897, 0.000722, 0.001248, 0.000668, 0.000585, 0.001174, 0.009094, 0.001413, 0.001731, 0.000866, 0.002782, 0.000702, 0.000377, 0.002367, 0.014317, 0.011161, 0.004303, 0.083316, 0.023953, 0.001174, 0.008709, 0.001468, 0.002953, 0.003776, 0.011470, 0.001361, 0.005929, -0.674727, 0.014580, 0.289399,
                          0.005134, 0.001612, 0.009549, 0.003780, 0.007256, 0.004386, 0.003892, 0.003650, 0.002703, 0.002099, 0.002294, 0.001057, 0.014394, 0.026127, 0.070423, 0.029953, 0.004552, 0.001303, 0.007416, 0.001937, 0.001483, 0.001460, 0.001658, 0.000831, 0.001244, 0.002063, 0.001943, 0.001452, 0.062751, 0.158996, 0.352069, 0.156779, 0.002326, 0.001100, 0.004314, 0.002122, 0.004262, 0.003845, 0.004216, 0.005074, 0.000813, 0.000716, 0.000864, 0.000555, 0.008558, 0.014517, 0.030840, 0.013473, 0.004420, 0.002884, 0.004626, 0.002053, 0.002568, 0.001777, 0.000441, 0.009051, 0.001523, 0.176854, 0.030303, -1.304410, 0.024070,
                          0.002251, 0.002443, 0.002279, 0.003053, 0.002903, 0.003136, 0.001739, 0.002095, 0.000480, 0.002546, 0.000343, 0.001933, 0.004630, 0.012038, 0.017844, 0.022139, 0.001431, 0.004439, 0.001592, 0.006699, 0.000860, 0.000519, 0.001155, 0.001232, 0.000324, 0.000767, 0.000553, 0.000376, 0.005670, 0.012297, 0.028200, 0.023794, 0.000966, 0.000694, 0.001388, 0.002479, 0.001660, 0.003343, 0.001790, 0.004686, 0.001831, 0.001327, 0.000932, 0.000642, 0.002593, 0.005737, 0.010501, 0.013995, 0.035159, 0.072228, 0.001266, 0.002683, 0.002048, 0.007672, 0.001252, 0.013918, 0.001164, 0.013919, 0.485610, 0.019433, -0.882676 };

    Db_matrix *cPi = new Db_matrix(char_as,"pi_char");
    Db_matrix *cQ  = new Db_matrix(char_as,char_as,"Q_char");


    for(int j=0;j<char_as;j++)
    {
        cPi->s(tmp_pi[j],j);
        for(int i=0;i<char_as;i++)
        {
            cQ->s(tmp_q[j*char_as+i],j,i);
        }
    }


    /************************************************************/


    if(parsimony_table != 0)
        delete parsimony_table;

    parsimony_table = new Int_matrix(char_fas,char_fas,"parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            if(i==j)
            {
                parsimony_table->s(i,i,j);
            }
            else
            {
                Codon_symbol *codon1 = &codon_symbols.at(i);
                Codon_symbol *codon2 = &codon_symbols.at(j);

                if(codon1->index == char_as) // 'NNN'
                {
                    parsimony_table->s(j,i,j);
                }

                else if(codon2->index == char_as) // 'NNN'
                {
                    parsimony_table->s(i,i,j);
                }

                else if(codon1->n_units == 1 && codon2->n_units == 1)
                {
                    int c1 = min(codon1->first_codon,codon2->first_codon);
                    int c2 = max(codon1->first_codon,codon2->first_codon);

                    int add = char_as-2;
                    int sum = char_as;
                    for(int loop=0;loop<c1;loop++)
                    {
                        sum += add;
                        add--;
                    }
                    sum += c2;

                    parsimony_table->s(sum,i,j);
                }

                else if(codon1->n_units == 1 && codon2->n_units == 2 &&
                        ( codon1->first_codon == codon2->first_codon ||
                            codon1->first_codon == codon2->second_codon ) )
                {
                    parsimony_table->s(codon1->first_codon,i,j);
                }

                else if(codon2->n_units == 1 && codon1->n_units == 2 &&
                        ( codon2->first_codon == codon1->first_codon ||
                            codon2->first_codon == codon1->second_codon ) )
                {
                    parsimony_table->s(codon2->first_codon,i,j);
                }

                else
                {
                    // search for max score in Q matrix

                    int m = codon1->first_codon;
                    int n = codon2->first_codon;

                    float maxQ = cQ->g(m,n);
                    int maxl1 = m;
                    int maxl2 = n;

                    if(codon2->n_units == 2)
                    {
                        m = codon1->first_codon;
                        n = codon2->second_codon;

                        if(cQ->g(m,n)>maxQ)
                        {
                            maxQ = cQ->g(m,n);
                            maxl1 = codon1->first_codon;
                            maxl2 = codon2->second_codon;
                        }
                    }

                    if(codon1->n_units == 2)
                    {
                        m = codon1->second_codon;
                        n = codon2->first_codon;

                        if(cQ->g(m,n)>maxQ)
                        {
                            maxQ = cQ->g(m,n);
                            maxl1 = codon1->second_codon;
                            maxl2 = codon2->first_codon;
                        }
                    }

                    if(codon1->n_units == 2 && codon2->n_units == 2)
                    {
                        m = codon1->second_codon;
                        n = codon2->second_codon;

                        if(cQ->g(m,n)>maxQ)
                        {
                            maxQ = cQ->g(m,n);
                            maxl1 = codon1->second_codon;
                            maxl2 = codon2->second_codon;
                        }
                    }


                    // find the corresponding ambiguity character
                    int c1 = min(maxl1,maxl2);
                    int c2 = max(maxl1,maxl2);

                    int add = char_as-2;
                    int sum = char_as;
                    for(int loop=0;loop<c1;loop++)
                    {
                        sum += add;
                        add--;
                    }
                    sum += c2;

                    parsimony_table->s(sum,i,j);
                }

            }
        } // for(j)
    } // for(i)

    delete cPi;
    delete cQ;

//    this->print_int_matrix(parsimony_table);


    /*********************/

    // for situation where the parent's state has been updated and child's state may need to be changed.
    // if child's state is included in parent's state, use parsimony_table to get the minimum overlap;
    // if child's state is not included, change has happened between child and parent and child is not updated.

    if(child_parsimony_table != 0)
        delete child_parsimony_table;

    child_parsimony_table = new Int_matrix(char_fas,char_fas,"child_parsimony_char");

    for(int i=0;i<char_fas;i++)
    {
        for(int j=0;j<char_fas;j++)
        {
            Codon_symbol *parent = &codon_symbols.at(i);
            Codon_symbol *child = &codon_symbols.at(j);

            // parent and child are idnetical
            if(i == j)
            {
                child_parsimony_table->s(j,i,j);
            }

            // one of them is NNN -> the other must have more or equal amount of information
            else if(parent->index == char_as) // 'NNN'
            {
                child_parsimony_table->s(j,i,j);
            }

            else if(child->index == char_as) // 'NNN'
            {
                child_parsimony_table->s(i,i,j);
            }

            // child is known: no change
            else if(child->n_units == 1)
            {
                child_parsimony_table->s(j,i,j);
            }

            // parent is known: only two comparisons
            else if(parent->n_units == 1 && parent->symbol != "NNN")
            {
                if( parent->first_codon == child->first_codon || parent->second_codon == child->first_codon)
                    child_parsimony_table->s(i,i,j);
                else
                    child_parsimony_table->s(j,i,j);
            }

//            else if(child->n_units == 1 && parent->n_units == 1)
//            {
//                child_parsimony_table->s(j,i,j);
//            }

//            else if(child->n_units == 1 && parent->n_units == 2)
//            {
//                if( parent->first_codon == child->first_codon || parent->second_codon == child->first_codon)
//                    child_parsimony_table->s(j,i,j);
//                else
//                    child_parsimony_table->s(-2,i,j);
//            }

//            // parent is known -> child must be the same

//            else if(parent->n_units == 1 && child->n_units == 2 && parent->symbol != "NNN")
//            {
//                if(parent->first_codon == child->first_codon || parent->first_codon == child->second_codon)
//                    child_parsimony_table->s(i,i,j);
//                else
//                    child_parsimony_table->s(-3,i,j);
//            }

            // both are known: four comparisons
            else
            {
                int c = -1;
                if(parent->first_codon == child->first_codon || parent->first_codon == child->second_codon)
                    c = parent->first_codon;
                else if(parent->second_codon == child->first_codon || parent->second_codon == child->second_codon)
                    c = parent->second_codon;

                child_parsimony_table->s(c,i,j);
            }
        }
    }

    if(mostcommon_table != 0)
        delete mostcommon_table;


    mostcommon_table = new Int_matrix(char_as,char_as,"mostcommon_char");
    for(int i=0;i<char_as;i++)
    {
        for(int j=0;j<char_as;j++)
        {
            mostcommon_table->s(j,i,j);
            if(tmp_pi[i]>tmp_pi[j])
                mostcommon_table->s(i,i,j);
        }
    }

    if(Settings::noise>5)
    {
        stringstream ss;

        ss<<"\nModel_factory::define_codon_alphabet(). Protein parsimony table.\n\n  ";
        for(int i=0;i<char_fas;i++)
            ss<<full_character_alphabet->at(i)<<" ";
        ss<<endl;

        for(int i=0;i<char_fas;i++)
        {
            ss<<full_character_alphabet->at(i)<<" ";
            for(int j=0;j<char_fas;j++)
            {
                ss<<full_character_alphabet->at(parsimony_table->g(i,j))<<" ";
            }
            ss<<endl;
        }
        ss<<endl;

        Log_output::write_out(ss.str(),6);
    }

    if(char_ambiguity != 0)
        delete char_ambiguity;

    char_ambiguity = new Db_matrix(1,1); // not needed?
    char_ambiguity->initialise(0);

}

/*******************************************/


/*
 * For debugging.
 */
void Model_factory::print_char_alphabet()
{
    stringstream ss;

    ss<<"\nModel_factory::print_char_alphabet()\n";
    for(unsigned int i=0;i<char_symbols.size();i++)
    {
        Char_symbol a = char_symbols.at(i);

        ss<<"Index "<<a.index<<"; symbol "<<a.symbol<<"; residues ("<<a.n_units<<"): ";
        for(int j=0;j<a.n_units;j++){
            int t = char_alphabet.find(a.residues.at(j));
            ss<<a.residues.at(j)<<"["<<t<<"] ";
        }
        ss<<endl;
    }

    Log_output::write_out(ss.str(),6);
}

void Model_factory::print_codon_alphabet()
{
    stringstream ss;

    ss<<"\nModel_factory::print_codon_alphabet()\n";
    for(unsigned int i=0;i<codon_symbols.size();i++)
    {
        Codon_symbol *a = &codon_symbols.at(i);

        ss<<"Index "<<a->index<<"; symbol "<<a->symbol<<"; residues ("<<a->n_units<<"): ";
        ss<<endl;
    }

    Log_output::write_out(ss.str(),6);
}

/*******************************************/


/*
 * Definition of DNA model: Q matrix and indel rate/prob.
 */

void Model_factory::dna_model(float *char_pi,Settings *st)
{
    float char_kappa = 2.0;
    float char_rho = 1.0;
    float ins_rate = 0.01;
    float del_rate = 0.01;
    float gap_ext = 0.8;
    float end_gap_ext = 0.95;
    float break_gap_ext = 0.99;

    if(st->is("pacbio"))
    {
        ins_rate = 1;
        del_rate = 1;
        gap_ext = 0.1;
    }

    if(st->is("char-kappa"))
        char_kappa =  st->get("char-kappa").as<float>();

    if(st->is("char-rho"))
        char_rho =  st->get("char-rho").as<float>();

    if(st->is("ins-rate"))
        ins_rate =  st->get("ins-rate").as<float>();

    if(st->is("del-rate"))
        del_rate =  st->get("del-rate").as<float>();

    if(st->is("indel-rate"))
        ins_rate =  del_rate =  st->get("indel-rate").as<float>();

    if(st->is("gap-extension"))
        gap_ext =  st->get("gap-extension").as<float>();

    if(st->is("end-gap-extension"))
        end_gap_ext =  st->get("end-gap-extension").as<float>();

    if(st->is("pair-read-gap-extension"))
        break_gap_ext =  st->get("pair-read-gap-extension").as<float>();


    this->dna_model(char_pi,char_kappa,char_rho,ins_rate,del_rate,gap_ext,end_gap_ext,break_gap_ext);
}

void Model_factory::dna_model(float* pi,float kappa, float rho,float ins_rate,float del_rate, float ext_prob, float end_ext_prob, float break_ext_prob)
{

    if (Settings::noise>4){
        stringstream ss;
        ss<<"DNA substitution model: base frequencies "
                <<pi[0]<<", "<<pi[1]<<", "<<pi[2]<<" and "<<pi[3]
                <<"; kappa "<<kappa<<" and rho "<<rho<<"; ins rate "<<ins_rate<<", del rate "<<del_rate<<" and ext. prob "<<ext_prob<<"."<<endl;
        Log_output::write_out(ss.str(),5);
    }

    char_ins_rate = ins_rate;
    char_del_rate = del_rate;
    char_ext_prob = ext_prob;
    char_end_ext_prob = end_ext_prob;
    char_break_ext_prob = break_ext_prob;


    if(Settings::noise > 4)
        print_char_alphabet();

    // pi and Q
    // Allocate space for P matrices
    //
    Db_matrix *charQ  = new Db_matrix(char_as,char_as,"Q_char");

    charPi = new Db_matrix(char_fas,"pi_char");


    for(int j=0;j<char_as;j++)
    {
        charPi->s(pi[j],j);
    }

    float ka = kappa/2.0;

    float piR = pi[0]+pi[2];
    float piY = pi[1]+pi[3];

    float beta = 1/(2*piR*piY*(1+ka));

    float alfaY = ( piR*piY*ka - pi[0]*pi[2] - pi[1]*pi[3]) /
                  ( (2+2*ka) * ( piY*pi[0]*pi[2] * rho+piR*pi[1]*pi[3] ) );

    float alfaR = rho*alfaY;

    /////////////////////////////
    /* filling up the Q matrix */
    double t;

    /*1st row*/
    t = beta*pi[1];
    charQ->s(t,0,1);

    t = alfaR*pi[2] / piR+beta*pi[2];
    charQ->s(t,0,2);

    t = beta*pi[3];
    charQ->s(t,0,3);

    t = 0-charQ->g(0,1)-charQ->g(0,2)-charQ->g(0,3);
    charQ->s(t,0,0);

    /*2nd row*/
    t = beta*pi[0];
    charQ->s(t,1,0);

    t = beta*pi[2];
    charQ->s(t,1,2);

    t = alfaY*pi[3]/piY+beta*pi[3];
    charQ->s(t,1,3);

    t = 0-charQ->g(1,0)-charQ->g(1,2)-charQ->g(1,3);
    charQ->s(t,1,1);

    /*3rd row*/
    t = alfaR*pi[0]/piR+beta*pi[0];
    charQ->s(t,2,0);

    t = beta*pi[1];
    charQ->s(t,2,1);

    t = beta*pi[3];
    charQ->s(t,2,3);

    t = 0-charQ->g(2,0)-charQ->g(2,1)-charQ->g(2,3);
    charQ->s(t,2,2);

    /*4th row*/
    t = beta*pi[0];
    charQ->s(t,3,0);

    t = alfaY*pi[1]/piY+beta*pi[1];
    charQ->s(t,3,1);

    t = beta*pi[2];
    charQ->s(t,3,2);

    t = 0-charQ->g(3,0)-charQ->g(3,1)-charQ->g(3,2);
    charQ->s(t,3,3);

    /* filling up the Q matrix */
    /////////////////////////////

    // Find eigenvalues and eigenvectors.
    //
    charU = new Db_matrix(char_as,char_as,"eigenvectors_1");
    charU->initialise();
    charV = new Db_matrix(char_as,char_as,"eigenvectors_2");
    charV->initialise();
    charRoot = new Db_matrix(char_as,"eigenvalues");
    charRoot->initialise();

    build_model(char_as,charPi,charQ,charU,charV,charRoot);

    if (Settings::noise > 4) {
        print_char_q_matrices(charQ);
    }

    delete charQ;

    character_alphabet = this->get_dna_character_alphabet();
    full_character_alphabet = this->get_dna_full_character_alphabet();

    ancestral_character_alphabet.clear();

    for(int i=0;i<(int)dna_full_char_alphabet.length();i++)
        ancestral_character_alphabet.push_back(std::string(1,dna_full_char_alphabet.at(i)));

}

/*******************************************/

void Model_factory::protein_model(Settings *st)
{
    float ins_rate = 0.05;
    if(st->is("ins-rate"))
        ins_rate =  st->get("ins-rate").as<float>();

    float del_rate = 0.05;
    if(st->is("del-rate"))
        del_rate =  st->get("del-rate").as<float>();

    if(st->is("indel-rate"))
        ins_rate =  del_rate =  st->get("indel-rate").as<float>();

    float gap_ext = 0.5;
    if(st->is("gap-extension"))
        gap_ext =  st->get("gap-extension").as<float>();

    float end_gap_ext = 0.75;
    if(st->is("end-gap-extension"))
        end_gap_ext =  st->get("end-gap-extension").as<float>();

    this->protein_model(ins_rate,del_rate,gap_ext,end_gap_ext);
}

void Model_factory::protein_model(float ins_rate,float del_rate, float ext_prob, float end_ext_prob)
{

    Log_output::write_out("Protein substitution model: WAG\n",5);

    char_ins_rate = ins_rate;
    char_del_rate = del_rate;
    char_ext_prob = ext_prob;
    char_end_ext_prob = end_ext_prob;
    char_break_ext_prob = 0.0;


    if(Settings::noise > 4)
        print_char_alphabet();

    // pi and Q
    // Allocate space for P matrices
    //
    double tmp_pi[20] = {0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466,
                         0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956};

    charPi = new Db_matrix(char_fas,"pi_char");

    for(int j=0;j<char_as;j++)
    {
        charPi->s(tmp_pi[j],j);
    }

    double tmp_q[400] = {-1.0644447077525, 0.024253680012, 0.0199296524112, 0.0421562148098, 0.019829882912, 0.0333710782038, 0.091898529865, 0.117944490096, 0.0077435982602, 0.00937017411, 0.034303854235, 0.056214349179, 0.0174255844392, 0.0080896843586, 0.065832507505, 0.234330242141, 0.129414648097, 0.0016275200247, 0.008491734537, 0.142217282556,
                         0.0477814374309, -0.9283507809291, 0.0248352939324, 0.0084029714104, 0.0101982061898, 0.11148814755, 0.0254969723473, 0.048674413647, 0.052213352795, 0.009062124214, 0.042903719239, 0.331941090612, 0.0133235035374, 0.0039473788809, 0.0310955230559, 0.085103118001, 0.0338262340451, 0.016744036728, 0.0134582713486, 0.0178549859644,
                         0.0441670615592, 0.027937434312, -1.38391347672512, 0.309721806842, 0.0051215097968, 0.056694964284, 0.0549932739622, 0.093704896008, 0.096657307877, 0.026861601976, 0.011338897352, 0.186830763486, 0.0038658446967, 0.00369569221099, 0.0089275113111, 0.276280123717, 0.123859441762, 0.00103458645453, 0.0383077812, 0.0139129779176,
                         0.0640178448442, 0.006477251488, 0.212232770148, -0.94297329251968, 0.00058492787022, 0.0226532677023, 0.358464938024, 0.0720614260512, 0.0227376245588, 0.001911353642, 0.0073109283823, 0.029764733853, 0.0020234831358, 0.00179593805976, 0.0194028221904, 0.074506504504, 0.0228715867982, 0.0018668150853, 0.0114891949562, 0.010799881226,
                         0.088970318416, 0.023225614652, 0.0103686978864, 0.00172817559999, -0.46436390434772, 0.00362939371299, 0.0012396736328, 0.0255311625132, 0.0060827096236, 0.00824576291, 0.033128997983, 0.00459221916954, 0.0076154533014, 0.015296664838, 0.0050066661924, 0.097857567114, 0.0312985388968, 0.010315697313, 0.0191832740086, 0.071047316584,
                         0.0787099366842, 0.133477006, 0.060339961416, 0.0351844479133, 0.00190795624962, -1.31468853172694, 0.317551411783, 0.0274774230936, 0.104910689643, 0.005521101322, 0.074957777201, 0.24159519414, 0.030136742202, 0.00384014619352, 0.0427139961732, 0.071524881773, 0.0523445036856, 0.0031035709083, 0.008032288082, 0.0213594972636,
                         0.137118971515, 0.019310611604, 0.0370254015012, 0.352205574616, 0.0004122601456, 0.200883241107, -1.17853838814971, 0.0472634621406, 0.0139264517825, 0.00617432607, 0.013298858967, 0.160308574698, 0.0061457688348, 0.00311812993141, 0.0312266801005, 0.0490058789081, 0.0501991141155, 0.0022522133463, 0.0069244312826, 0.0417384374836,
                         0.122727478488, 0.02570888938, 0.043997465064, 0.0493773258384, 0.0059212002572, 0.0121221828612, 0.0329610245313, -0.4741383934945, 0.006093410533, 0.0014757945466, 0.0052849306733, 0.0231712797588, 0.00339542007, 0.0019189431989, 0.011146518267, 0.093280508578, 0.0137786810791, 0.0048478037397, 0.0036545482168, 0.0132749884132,
                         0.0274570594166, 0.0939747598, 0.154649002326, 0.0530905054876, 0.0048071015816, 0.157714501491, 0.0330950244725, 0.020763831438, -0.9455247927498, 0.00669751654, 0.043058119558, 0.0552322503552, 0.0078818406807, 0.0261095183349, 0.0318601786938, 0.0514549945251, 0.0288777379989, 0.0037772913771, 0.136632497248, 0.0083910614248,
                         0.0167482050465, 0.008221840588, 0.0216647526984, 0.0022496876087, 0.003284932553, 0.0041839549677, 0.0073964135655, 0.00253502563518, 0.003376161347, -1.17498318692916, 0.27336615273, 0.0200868455952, 0.083031965142, 0.040717445093, 0.00457305166728, 0.022206797976, 0.088966278632, 0.0030567591897, 0.014821160614, 0.55449575628,
                         0.0344705408285, 0.021883589212, 0.0051413506032, 0.00483769259197, 0.0074197365386, 0.0319346789409, 0.0089563400907, 0.00510364337166, 0.0122025059606, 0.15368423202, -0.69175870029473, 0.015975776073, 0.094666495854, 0.081290001923, 0.0190303105564, 0.0239655313281, 0.0199280900994, 0.0095710687431, 0.0140609310556, 0.127636184504,
                         0.0785078337935, 0.23531264024, 0.117737663694, 0.0273733764605, 0.00142943173442, 0.14305227669, 0.150049162927, 0.0310993759044, 0.0217544113216, 0.015694841712, 0.022203558995, -1.07152398375532, 0.0182209045452, 0.0034141362684, 0.0254852873376, 0.067232846627, 0.084623394646, 0.0019781331795, 0.0047007809888, 0.0216539266904,
                         0.0774016821384, 0.030039999464, 0.0077483399574, 0.0059186573054, 0.0075393483596, 0.056754463806, 0.0182957528036, 0.01449413838, 0.0098736900133, 0.20634205636, 0.41846021018, 0.0579518322936, -1.2597234186325, 0.045758173097, 0.0078405461599, 0.0343352383995, 0.092502574724, 0.0074188949454, 0.0151127724254, 0.14593504782,
                         0.0182346531826, 0.004516408092, 0.00375891879174, 0.00266574034104, 0.007684890556, 0.00366990113448, 0.00471054498671, 0.0041568456258, 0.0165979167123, 0.05134827302, 0.18234669053, 0.0055103727096, 0.023220499701, -0.67999937105987, 0.0073881779164, 0.0379519766649, 0.0104882661681, 0.022005248076, 0.227669563576, 0.0460744832752,
                         0.124618565545, 0.029878490308, 0.0076255992414, 0.0241862096784, 0.0021123505512, 0.0342809801532, 0.0396167807095, 0.020277640926, 0.0170090221974, 0.0048431492208, 0.035849495396, 0.0345434792256, 0.0033413780883, 0.0062045996636, -0.5770185229931, 0.112151837712, 0.0485285253768, 0.0020054663895, 0.0076208498132, 0.0223241027972,
                         0.292004459041, 0.05383008268, 0.155350266162, 0.061138656376, 0.027178817748, 0.037788440247, 0.0409279829071, 0.111708930276, 0.0180832908897, 0.01548197904, 0.029719604451, 0.059989719918, 0.0096324810435, 0.0209811655989, 0.073828693968, -1.326554610767, 0.267114820854, 0.0075345000378, 0.0277605484806, 0.0165001710484,
                         0.183747304969, 0.024378648436, 0.079353827364, 0.0213842684566, 0.0099045924752, 0.0315100653768, 0.0477688308585, 0.0188010037494, 0.0115635053091, 0.07067118256, 0.028157755998, 0.086032427628, 0.029568433524, 0.0066065589057, 0.0363992375304, 0.304350756558, -1.1004826896859, 0.0015948784176, 0.0102700127816, 0.098419398788,
                         0.0098004742107, 0.05117989024, 0.00281118065298, 0.0074025714917, 0.013845044146, 0.0079236101097, 0.0090895272073, 0.0280544413194, 0.0064149020097, 0.010298201078, 0.057355623581, 0.008529242643, 0.0100576594062, 0.058786971516, 0.0063796049555, 0.0364094439818, 0.0067641119728, -0.44467569893618, 0.087670143938, 0.0259030544764,
                         0.0208543675065, 0.016776769076, 0.0424510884, 0.0185802165661, 0.0105002187974, 0.008363355651, 0.0113971362467, 0.0086252194872, 0.094633174672, 0.02036395922, 0.034364459162, 0.0082661793504, 0.0083556782799, 0.248050243532, 0.0098869347026, 0.0547101006747, 0.0177637255796, 0.035754572001, -0.6920103710931, 0.022312972188,
                         0.173776433679, 0.011074304228, 0.0076711383924, 0.0086899653085, 0.019349118692, 0.0110654786961, 0.0341810742559, 0.0155886497946, 0.0028916398054, 0.3790671258, 0.15520551106, 0.0189456434124, 0.040145332815, 0.0249765843548, 0.0144102052697, 0.0161795265281, 0.084699660521, 0.0052561618971, 0.011101848966, -1.034275403476 };

    Db_matrix *charQ  = new Db_matrix(char_as,char_as,"Q_char");

    for(int j=0;j<char_as;j++)
    {
        for(int i=0;i<char_as;i++)
        {
            charQ->s(tmp_q[j*char_as+i],j,i);
        }
    }

    // Find eigenvalues and eigenvectors.
    //
    charU = new Db_matrix(char_as,char_as,"eigenvectors_1");
    charU->initialise();
    charV = new Db_matrix(char_as,char_as,"eigenvectors_2");
    charV->initialise();
    charRoot = new Db_matrix(char_as,"eigenvalues");
    charRoot->initialise();

    build_model(char_as,charPi,charQ,charU,charV,charRoot);

    if (Settings::noise > 4) {
        print_char_q_matrices(charQ);
    }

    delete charQ;

    character_alphabet = this->get_protein_character_alphabet();
    full_character_alphabet = this->get_protein_full_character_alphabet();

    ancestral_character_alphabet.clear();

    for(int i=0;i<(int)protein_char_alphabet.length();i++)
        ancestral_character_alphabet.push_back(std::string(1,protein_char_alphabet.at(i)));

    ancestral_character_alphabet.push_back("X");

    for(int i=0;i<(int)protein_char_alphabet.length()-1;i++)
        for(int j=i+1;j<(int)protein_char_alphabet.length();j++)
            if(charPi->g(i) > charPi->g(j))
                ancestral_character_alphabet.push_back(std::string(1,protein_char_alphabet.at(i)));
            else
                ancestral_character_alphabet.push_back(std::string(1,protein_char_alphabet.at(j)));

}

/*******************************************/

void Model_factory::codon_model(Settings *st)
{
    float ins_rate = 0.01;
    if(st->is("ins-rate"))
        ins_rate =  st->get("ins-rate").as<float>();

    float del_rate = 0.01;
    if(st->is("del-rate"))
        del_rate =  st->get("del-rate").as<float>();

    if(st->is("indel-rate"))
        ins_rate =  del_rate =  st->get("indel-rate").as<float>();

    float gap_ext = 0.5;
    if(st->is("gap-extension"))
        gap_ext =  st->get("gap-extension").as<float>();

    float end_gap_ext = 0.75;
    if(st->is("end-gap-extension"))
        end_gap_ext =  st->get("end-gap-extension").as<float>();

    this->codon_model(ins_rate,del_rate,gap_ext,end_gap_ext);
}


void Model_factory::codon_model(float ins_rate,float del_rate, float ext_prob, float end_ext_prob)
{

    Log_output::write_out("Codon substitution model: Kosiol&Goldman\n",5);

    char_ins_rate = ins_rate;
    char_del_rate = del_rate;
    char_ext_prob = ext_prob;
    char_end_ext_prob = end_ext_prob;
    char_break_ext_prob = 0.0;

    if(Settings::noise > 4)
        print_char_alphabet();

      double tmp_pi[61] = {0.024709, 0.026990, 0.037138, 0.018760, 0.013945, 0.024015, 0.009522, 0.014230, 0.010544, 0.015952, 0.008214, 0.008979, 0.007064, 0.029436, 0.020724, 0.019431, 0.015416, 0.013015, 0.026398, 0.008551, 0.011465, 0.012756, 0.006648, 0.009507, 0.005177, 0.012128, 0.007718, 0.007807, 0.005328, 0.017807, 0.031660, 0.011119, 0.028694, 0.027083, 0.035815, 0.023890, 0.013035, 0.027437, 0.008784, 0.019503, 0.015756, 0.021140, 0.009781, 0.015576, 0.007598, 0.019661, 0.026815, 0.015815, 0.022021, 0.013288, 0.009289, 0.015669, 0.007124, 0.012150, 0.019672, 0.015708, 0.014291, 0.007146, 0.026531, 0.012766, 0.015811};

    charPi = new Db_matrix(char_fas,"pi_char");

    for(int j=0;j<char_as;j++)
    {
        charPi->s(tmp_pi[j],j);
    }

    double tmp_q[3721] = {-0.996852, 0.019663, 0.434175, 0.025324, 0.017695, 0.012761, 0.005660, 0.011280, 0.061068, 0.009776, 0.021005, 0.009536, 0.005768, 0.002711, 0.011898, 0.003740, 0.042628, 0.005906, 0.026732, 0.007654, 0.007352, 0.002504, 0.002391, 0.003500, 0.012882, 0.015172, 0.009349, 0.014707, 0.001452, 0.001182, 0.001249, 0.003244, 0.061214, 0.009708, 0.018164, 0.012120, 0.009128, 0.004818, 0.003400, 0.009158, 0.004921, 0.003252, 0.000831, 0.002467, 0.003219, 0.004301, 0.002816, 0.003567, 0.001069, 0.002186, 0.008633, 0.004091, 0.002042, 0.006228, 0.000140, 0.001428, 0.000110, 0.002805, 0.000981, 0.002652, 0.001441,
                          0.018001, -0.833627, 0.030232, 0.299934, 0.007529, 0.029937, 0.005116, 0.010141, 0.003332, 0.059604, 0.002731, 0.016622, 0.001258, 0.003052, 0.006465, 0.001080, 0.008318, 0.033977, 0.017907, 0.012141, 0.002650, 0.004096, 0.001733, 0.002844, 0.002169, 0.005276, 0.001991, 0.003428, 0.000548, 0.001018, 0.001301, 0.000545, 0.009888, 0.070669, 0.013486, 0.023938, 0.002984, 0.013093, 0.003870, 0.004338, 0.008120, 0.023487, 0.007724, 0.006311, 0.001142, 0.002033, 0.001953, 0.001510, 0.007353, 0.001702, 0.005085, 0.012603, 0.004073, 0.005464, 0.001329, 0.000791, 0.000847, 0.000644, 0.002019, 0.000763, 0.001431,
                          0.288865, 0.021971, -0.843432, 0.014913, 0.007281, 0.018314, 0.011510, 0.007363, 0.023323, 0.011365, 0.053669, 0.004764, 0.001533, 0.002715, 0.014706, 0.001548, 0.016566, 0.007204, 0.058416, 0.005632, 0.003355, 0.003986, 0.003516, 0.002606, 0.007632, 0.031929, 0.024684, 0.011513, 0.001164, 0.002390, 0.004846, 0.001680, 0.018442, 0.013350, 0.042174, 0.006639, 0.004700, 0.008692, 0.007474, 0.007676, 0.004622, 0.006910, 0.002998, 0.001889, 0.002176, 0.004544, 0.005439, 0.002552, 0.002017, 0.002332, 0.003315, 0.006478, 0.006181, 0.004989, 0.000223, 0.000560, 0.000073, 0.000814, 0.000961, 0.003282, 0.000970,
                          0.033355, 0.431526, 0.029523, -0.987099, 0.014903, 0.010645, 0.007384, 0.023703, 0.006355, 0.020501, 0.004741, 0.036958, 0.001072, 0.001118, 0.006785, 0.003785, 0.020131, 0.010786, 0.014616, 0.024978, 0.004115, 0.002128, 0.001443, 0.003259, 0.002156, 0.001891, 0.002442, 0.003918, 0.000564, 0.000919, 0.001150, 0.002896, 0.020626, 0.031842, 0.011045, 0.072138, 0.004404, 0.004479, 0.003079, 0.007530, 0.012646, 0.007896, 0.006704, 0.014683, 0.001177, 0.001220, 0.002259, 0.004108, 0.003652, 0.006112, 0.009177, 0.006368, 0.004857, 0.011702, 0.000900, 0.000509, 0.001259, 0.000831, 0.001006, 0.002572, 0.002573,
                          0.031354, 0.014573, 0.019391, 0.020049, -1.243500, 0.259288, 0.147445, 0.260513, 0.015627, 0.024450, 0.004176, 0.020780, 0.021256, 0.005638, 0.021888, 0.008172, 0.016821, 0.003028, 0.012152, 0.006987, 0.007839, 0.003930, 0.004451, 0.003892, 0.002448, 0.002018, 0.003008, 0.002370, 0.002550, 0.002950, 0.004344, 0.003024, 0.017739, 0.005225, 0.010136, 0.006291, 0.021145, 0.009424, 0.008924, 0.012785, 0.003689, 0.002751, 0.001829, 0.002262, 0.011715, 0.010648, 0.016584, 0.011067, 0.002333, 0.003536, 0.048425, 0.021486, 0.017250, 0.020995, 0.000984, 0.000529, 0.001355, 0.003769, 0.004284, 0.006642, 0.003291,
                          0.013129, 0.033646, 0.028321, 0.008316, 0.150556, -1.097510, 0.139816, 0.197266, 0.002217, 0.059113, 0.002412, 0.014527, 0.003885, 0.034135, 0.020238, 0.004692, 0.007230, 0.009949, 0.019241, 0.003572, 0.001964, 0.007404, 0.002603, 0.002921, 0.001748, 0.009255, 0.003360, 0.001426, 0.000805, 0.009008, 0.006028, 0.001764, 0.009818, 0.013657, 0.014062, 0.004799, 0.006814, 0.044137, 0.007396, 0.009147, 0.001874, 0.005363, 0.002093, 0.001265, 0.005544, 0.032109, 0.015770, 0.007693, 0.006001, 0.001725, 0.008444, 0.051019, 0.016328, 0.014919, 0.002173, 0.001355, 0.000560, 0.001103, 0.007393, 0.002331, 0.002065,
                          0.014687, 0.014502, 0.044894, 0.014547, 0.215929, 0.352630, -1.406200, 0.189339, 0.004280, 0.030649, 0.007156, 0.016698, 0.011924, 0.013213, 0.028368, 0.008267, 0.009388, 0.005630, 0.030042, 0.005121, 0.003732, 0.005286, 0.006377, 0.003060, 0.002635, 0.004339, 0.005676, 0.002939, 0.001577, 0.006201, 0.010942, 0.002943, 0.010931, 0.012754, 0.021891, 0.008278, 0.009608, 0.026715, 0.021540, 0.010194, 0.002835, 0.004400, 0.003624, 0.001051, 0.006111, 0.017671, 0.037522, 0.009466, 0.006702, 0.003924, 0.014147, 0.026283, 0.036125, 0.017230, 0.003900, 0.002086, 0.001308, 0.001115, 0.007715, 0.005218, 0.002888,
                          0.019586, 0.019234, 0.019218, 0.031249, 0.255292, 0.332923, 0.126698, -1.266280, 0.002834, 0.026528, 0.003170, 0.034494, 0.004782, 0.010430, 0.010621, 0.022796, 0.012602, 0.005524, 0.009239, 0.007390, 0.002554, 0.002426, 0.001468, 0.006049, 0.001136, 0.002621, 0.001335, 0.002962, 0.001345, 0.002380, 0.002930, 0.002783, 0.016675, 0.005733, 0.010980, 0.012087, 0.009933, 0.016212, 0.005285, 0.027461, 0.003338, 0.002825, 0.002368, 0.003573, 0.007188, 0.008487, 0.010058, 0.023212, 0.002478, 0.002207, 0.023887, 0.020905, 0.010822, 0.041072, 0.000922, 0.001716, 0.002281, 0.001789, 0.002582, 0.003275, 0.002328,
                          0.143109, 0.008529, 0.082151, 0.011307, 0.020667, 0.005050, 0.003865, 0.003825, -1.172030, 0.013413, 0.230080, 0.008486, 0.004219, 0.000918, 0.008596, 0.001573, 0.019807, 0.004337, 0.017706, 0.006355, 0.004002, 0.001226, 0.001383, 0.002244, 0.119639, 0.086404, 0.077589, 0.187931, 0.002005, 0.001361, 0.002177, 0.002336, 0.011718, 0.002232, 0.005356, 0.004616, 0.004753, 0.003025, 0.001352, 0.004934, 0.009323, 0.002279, 0.003855, 0.001965, 0.003069, 0.001752, 0.002632, 0.002507, 0.001718, 0.002002, 0.005187, 0.001598, 0.001824, 0.002826, 0.000740, 0.001169, 0.000452, 0.002036, 0.000829, 0.003273, 0.000720,
                          0.015143, 0.100850, 0.026460, 0.024110, 0.021374, 0.088994, 0.018295, 0.023664, 0.008865, -1.134840, 0.011222, 0.184047, 0.001266, 0.005081, 0.008771, 0.002019, 0.010616, 0.015081, 0.025476, 0.005128, 0.005533, 0.009222, 0.004831, 0.005342, 0.004001, 0.018797, 0.006563, 0.001963, 0.001404, 0.003401, 0.005569, 0.000883, 0.013753, 0.044440, 0.024604, 0.010764, 0.006736, 0.035682, 0.008081, 0.007624, 0.010821, 0.047683, 0.014155, 0.010143, 0.001711, 0.006961, 0.005633, 0.002857, 0.007246, 0.002259, 0.023880, 0.083006, 0.042204, 0.035529, 0.012829, 0.001691, 0.002213, 0.000597, 0.003560, 0.001680, 0.002524,
                          0.063186, 0.008973, 0.242659, 0.010827, 0.007089, 0.007053, 0.008295, 0.005492, 0.295346, 0.021795, -1.428430, 0.013076, 0.002855, 0.002348, 0.026156, 0.001463, 0.013660, 0.008049, 0.052049, 0.005814, 0.002194, 0.005230, 0.002925, 0.002948, 0.064463, 0.136976, 0.196339, 0.077079, 0.000977, 0.002183, 0.008981, 0.002248, 0.004789, 0.005236, 0.016869, 0.002400, 0.004462, 0.006344, 0.005525, 0.004035, 0.004823, 0.004205, 0.015277, 0.002281, 0.001683, 0.001982, 0.010400, 0.001254, 0.003289, 0.003215, 0.002211, 0.004172, 0.003565, 0.002322, 0.001597, 0.009974, 0.000339, 0.000816, 0.002412, 0.003565, 0.000659,
                          0.026242, 0.049966, 0.019705, 0.077218, 0.032273, 0.038856, 0.017708, 0.054666, 0.009965, 0.326984, 0.011962, -1.277050, 0.002464, 0.002193, 0.006083, 0.007151, 0.017106, 0.008477, 0.026180, 0.010608, 0.004151, 0.002994, 0.001776, 0.005592, 0.003719, 0.006082, 0.002849, 0.004285, 0.001048, 0.001740, 0.003428, 0.002061, 0.022995, 0.022963, 0.018742, 0.030636, 0.017795, 0.016301, 0.008316, 0.021800, 0.015141, 0.022784, 0.014023, 0.019830, 0.001705, 0.002613, 0.003840, 0.002359, 0.002096, 0.006557, 0.051613, 0.050652, 0.044803, 0.075464, 0.003223, 0.001183, 0.006296, 0.000676, 0.002201, 0.001503, 0.003405,
                          0.020177, 0.004806, 0.008062, 0.002846, 0.041962, 0.013209, 0.016074, 0.009634, 0.006297, 0.002858, 0.003319, 0.003132, -1.366290, 0.307296, 0.071081, 0.317317, 0.004254, 0.002661, 0.004990, 0.002242, 0.001804, 0.001521, 0.002119, 0.002147, 0.002004, 0.001584, 0.001030, 0.001793, 0.012765, 0.020674, 0.042996, 0.013209, 0.010283, 0.001041, 0.004226, 0.000884, 0.008055, 0.005836, 0.003125, 0.005231, 0.000923, 0.001197, 0.000587, 0.000291, 0.079194, 0.042651, 0.087968, 0.059847, 0.003415, 0.004866, 0.003778, 0.002968, 0.001696, 0.002722, 0.001694, 0.001694, 0.000923, 0.037165, 0.011791, 0.026013, 0.010364,
                          0.002275, 0.002798, 0.003425, 0.000712, 0.002671, 0.027849, 0.004274, 0.005042, 0.000329, 0.002754, 0.000655, 0.000669, 0.073742, -0.887525, 0.036537, 0.277085, 0.000588, 0.001591, 0.002204, 0.000938, 0.000462, 0.000916, 0.000432, 0.000442, 0.000310, 0.001433, 0.000481, 0.000205, 0.005975, 0.051463, 0.049826, 0.011136, 0.000537, 0.001757, 0.001098, 0.000491, 0.001488, 0.008330, 0.001383, 0.001952, 0.000316, 0.000749, 0.000322, 0.000261, 0.011818, 0.116456, 0.077458, 0.032484, 0.003766, 0.001830, 0.000435, 0.001592, 0.000552, 0.000764, 0.003157, 0.001140, 0.001163, 0.005466, 0.023741, 0.011331, 0.006466,
                          0.014186, 0.008420, 0.026354, 0.006142, 0.014728, 0.023453, 0.013034, 0.007293, 0.004374, 0.006752, 0.010367, 0.002635, 0.024228, 0.051896, -0.819425, 0.034272, 0.009840, 0.004428, 0.033136, 0.005408, 0.001793, 0.002788, 0.002853, 0.001316, 0.001677, 0.003878, 0.003750, 0.001634, 0.015012, 0.042873, 0.100407, 0.022120, 0.005320, 0.001841, 0.010676, 0.001606, 0.007612, 0.017861, 0.009238, 0.012917, 0.003152, 0.004325, 0.003224, 0.000943, 0.009435, 0.016533, 0.059619, 0.012872, 0.009247, 0.003839, 0.004118, 0.006914, 0.003795, 0.002646, 0.005194, 0.007425, 0.001392, 0.017283, 0.022359, 0.043380, 0.013614,
                          0.004756, 0.001500, 0.002959, 0.003654, 0.005865, 0.005799, 0.004051, 0.016694, 0.000854, 0.001658, 0.000618, 0.003304, 0.115351, 0.419742, 0.036551, -1.067930, 0.002964, 0.001118, 0.001734, 0.002224, 0.000596, 0.000262, 0.000268, 0.001307, 0.000385, 0.000430, 0.000447, 0.000391, 0.009769, 0.015524, 0.031405, 0.039917, 0.002721, 0.000739, 0.001379, 0.001514, 0.004182, 0.003910, 0.001435, 0.006983, 0.000798, 0.000572, 0.000571, 0.000505, 0.021872, 0.045457, 0.058521, 0.112609, 0.001667, 0.002998, 0.000872, 0.000815, 0.000619, 0.002007, 0.001237, 0.000746, 0.001347, 0.014550, 0.007522, 0.019678, 0.018014,
                          0.068322, 0.014562, 0.039907, 0.024496, 0.015215, 0.011262, 0.005798, 0.011632, 0.013547, 0.010984, 0.007278, 0.009963, 0.001949, 0.001123, 0.013227, 0.003736, -1.042850, 0.017118, 0.362704, 0.020271, 0.012857, 0.002382, 0.004217, 0.004555, 0.018466, 0.010309, 0.008920, 0.008461, 0.003882, 0.002895, 0.006793, 0.007226, 0.096006, 0.010655, 0.041474, 0.019065, 0.012615, 0.006118, 0.003430, 0.011218, 0.007635, 0.002822, 0.002564, 0.003833, 0.003724, 0.003495, 0.006193, 0.004666, 0.003032, 0.005180, 0.010133, 0.010580, 0.004960, 0.012746, 0.000311, 0.004037, 0.000749, 0.005155, 0.001154, 0.003769, 0.001467,
                          0.011212, 0.070462, 0.020558, 0.015546, 0.003244, 0.018358, 0.004119, 0.006040, 0.003514, 0.018484, 0.005080, 0.005848, 0.001444, 0.003598, 0.007051, 0.001670, 0.020276, -0.982175, 0.086164, 0.296699, 0.002026, 0.007498, 0.005144, 0.002620, 0.005918, 0.041083, 0.014976, 0.006844, 0.001322, 0.012436, 0.006408, 0.002790, 0.009192, 0.032040, 0.018704, 0.010491, 0.002675, 0.007635, 0.002781, 0.002644, 0.002334, 0.007338, 0.003462, 0.001320, 0.001756, 0.007376, 0.005956, 0.001372, 0.069267, 0.022290, 0.003399, 0.013588, 0.005789, 0.006074, 0.001508, 0.005598, 0.001004, 0.000925, 0.020558, 0.001278, 0.005392,
                          0.025021, 0.018309, 0.082182, 0.010387, 0.006419, 0.017504, 0.010836, 0.004980, 0.007072, 0.015395, 0.016195, 0.008904, 0.001335, 0.002457, 0.026013, 0.001277, 0.211815, 0.042481, -0.912145, 0.019229, 0.005688, 0.005581, 0.005953, 0.005359, 0.007026, 0.023328, 0.028657, 0.005375, 0.002017, 0.007659, 0.022625, 0.003488, 0.025534, 0.012897, 0.104029, 0.008146, 0.004521, 0.013093, 0.009729, 0.005362, 0.003817, 0.006040, 0.006087, 0.000923, 0.001709, 0.002517, 0.008074, 0.001461, 0.005150, 0.002767, 0.004550, 0.009790, 0.008290, 0.004397, 0.000303, 0.003554, 0.000442, 0.001441, 0.002415, 0.003586, 0.000954,
                          0.022117, 0.038319, 0.024460, 0.054797, 0.011395, 0.010031, 0.005702, 0.012297, 0.007836, 0.009566, 0.005585, 0.011138, 0.001852, 0.003228, 0.013105, 0.005053, 0.036545, 0.451575, 0.059362, -1.154830, 0.005724, 0.007304, 0.002793, 0.011668, 0.007607, 0.013909, 0.012681, 0.014874, 0.001300, 0.006561, 0.004202, 0.006791, 0.012162, 0.014662, 0.014606, 0.017232, 0.004931, 0.006277, 0.006123, 0.004896, 0.008044, 0.005020, 0.005979, 0.002388, 0.003028, 0.003561, 0.008095, 0.002696, 0.030126, 0.057219, 0.008048, 0.007923, 0.003978, 0.013302, 0.001946, 0.006667, 0.000645, 0.003330, 0.007293, 0.002892, 0.012386,
                          0.015844, 0.006239, 0.010867, 0.006734, 0.009535, 0.004114, 0.003100, 0.003170, 0.003680, 0.007699, 0.001572, 0.003251, 0.001111, 0.001186, 0.003241, 0.001010, 0.017288, 0.002300, 0.013097, 0.004270, -0.959559, 0.223948, 0.145555, 0.273040, 0.003185, 0.001680, 0.003125, 0.003955, 0.002936, 0.003838, 0.004289, 0.001923, 0.018380, 0.008443, 0.014106, 0.011581, 0.018027, 0.005523, 0.006333, 0.009839, 0.003453, 0.003361, 0.002370, 0.001994, 0.001710, 0.001503, 0.002242, 0.003604, 0.001639, 0.001006, 0.026870, 0.006940, 0.005643, 0.011351, 0.000348, 0.001150, 0.000191, 0.001656, 0.000677, 0.001652, 0.001186,
                          0.004851, 0.008666, 0.011606, 0.003130, 0.004296, 0.013940, 0.003946, 0.002706, 0.001014, 0.011533, 0.003368, 0.002108, 0.000842, 0.002114, 0.004529, 0.000399, 0.002879, 0.007651, 0.011550, 0.004897, 0.201287, -0.935790, 0.169605, 0.242298, 0.001797, 0.008043, 0.004252, 0.002579, 0.001284, 0.014390, 0.006850, 0.002471, 0.005843, 0.011779, 0.014689, 0.003705, 0.005502, 0.032614, 0.005974, 0.005071, 0.002405, 0.007984, 0.002819, 0.001219, 0.001962, 0.006785, 0.004305, 0.002786, 0.002103, 0.001358, 0.005064, 0.029821, 0.006204, 0.006575, 0.000827, 0.000599, 0.000124, 0.000385, 0.004301, 0.001462, 0.000644,
                          0.008887, 0.007035, 0.019642, 0.004072, 0.009336, 0.009404, 0.009134, 0.003141, 0.002194, 0.011591, 0.003614, 0.002398, 0.002252, 0.001915, 0.008893, 0.000783, 0.009780, 0.010071, 0.023640, 0.003593, 0.251028, 0.325436, -1.182370, 0.186111, 0.004204, 0.006686, 0.005666, 0.003038, 0.002577, 0.004869, 0.011430, 0.002942, 0.014556, 0.014253, 0.023689, 0.006993, 0.013193, 0.017819, 0.025836, 0.008778, 0.003664, 0.010686, 0.006029, 0.001110, 0.002483, 0.005086, 0.008847, 0.002354, 0.001142, 0.001789, 0.009471, 0.012101, 0.015888, 0.008013, 0.001294, 0.002566, 0.000513, 0.001184, 0.001734, 0.003184, 0.002748,
                          0.009096, 0.008075, 0.010179, 0.006431, 0.005709, 0.007379, 0.003065, 0.009054, 0.002489, 0.008963, 0.002547, 0.005281, 0.001595, 0.001369, 0.002868, 0.002672, 0.007386, 0.003586, 0.014879, 0.010495, 0.329258, 0.325081, 0.130133, -1.102160, 0.002969, 0.002771, 0.001888, 0.005911, 0.002293, 0.002767, 0.003912, 0.009884, 0.009055, 0.009660, 0.007452, 0.012225, 0.009565, 0.008649, 0.003323, 0.021290, 0.003185, 0.004371, 0.003693, 0.002225, 0.000986, 0.003033, 0.003503, 0.003594, 0.002146, 0.003193, 0.007029, 0.006857, 0.007335, 0.022330, 0.000201, 0.001501, 0.000337, 0.001061, 0.001215, 0.001115, 0.002049,
                          0.061482, 0.011309, 0.054753, 0.007813, 0.006592, 0.008107, 0.004847, 0.003122, 0.243662, 0.012327, 0.102276, 0.006450, 0.002734, 0.001762, 0.006711, 0.001444, 0.054989, 0.014878, 0.035828, 0.012564, 0.007053, 0.004428, 0.005399, 0.005453, -1.458750, 0.251947, 0.223508, 0.179170, 0.004224, 0.002538, 0.006210, 0.003035, 0.006631, 0.005195, 0.011259, 0.003910, 0.005948, 0.006714, 0.002751, 0.003063, 0.005696, 0.005783, 0.003366, 0.003023, 0.003245, 0.002399, 0.003680, 0.001829, 0.002609, 0.001321, 0.007931, 0.004985, 0.004540, 0.005820, 0.000648, 0.001916, 0.000737, 0.001275, 0.001801, 0.003068, 0.000989,
                          0.030911, 0.011743, 0.097776, 0.002925, 0.002320, 0.018328, 0.003406, 0.003076, 0.075119, 0.024724, 0.092770, 0.004503, 0.000922, 0.003479, 0.006626, 0.000689, 0.013104, 0.044088, 0.050777, 0.009807, 0.001589, 0.008459, 0.003665, 0.002172, 0.107549, -1.141450, 0.182320, 0.204540, 0.001315, 0.010476, 0.007677, 0.002957, 0.004510, 0.005214, 0.012346, 0.002033, 0.001404, 0.010192, 0.002825, 0.002443, 0.001313, 0.005931, 0.001950, 0.001221, 0.001125, 0.006989, 0.003346, 0.001499, 0.010181, 0.001788, 0.001883, 0.011372, 0.004751, 0.002943, 0.002711, 0.004452, 0.000215, 0.000722, 0.003107, 0.002172, 0.001000,
                          0.029932, 0.006962, 0.118786, 0.005935, 0.005434, 0.010457, 0.007003, 0.002461, 0.106003, 0.013565, 0.208965, 0.003315, 0.000942, 0.001835, 0.010070, 0.001125, 0.017819, 0.025256, 0.098023, 0.014051, 0.004642, 0.007029, 0.004881, 0.002326, 0.149932, 0.286510, -1.403230, 0.121475, 0.001491, 0.003138, 0.015265, 0.002298, 0.003852, 0.003249, 0.020813, 0.001939, 0.002406, 0.007198, 0.003845, 0.004052, 0.002357, 0.005332, 0.006612, 0.001242, 0.001155, 0.001502, 0.009330, 0.001263, 0.002561, 0.001966, 0.002016, 0.003848, 0.005643, 0.001572, 0.003441, 0.009626, 0.001129, 0.001244, 0.002761, 0.003215, 0.001132,
                          0.046550, 0.011851, 0.054770, 0.009416, 0.004233, 0.004387, 0.003585, 0.005400, 0.253828, 0.004012, 0.081101, 0.004928, 0.001622, 0.000774, 0.004338, 0.000973, 0.016709, 0.011410, 0.018176, 0.016293, 0.005808, 0.004214, 0.002587, 0.007199, 0.118820, 0.317765, 0.120091, -1.214070, 0.001018, 0.001176, 0.003451, 0.004755, 0.005100, 0.002322, 0.006418, 0.003142, 0.002528, 0.003315, 0.001071, 0.003160, 0.001590, 0.002880, 0.001780, 0.001881, 0.001542, 0.002053, 0.002506, 0.001905, 0.004527, 0.003029, 0.002849, 0.002085, 0.002827, 0.005805, 0.000532, 0.000421, 0.001925, 0.001248, 0.001258, 0.002375, 0.000762,
                          0.006733, 0.002776, 0.008117, 0.001988, 0.006675, 0.003629, 0.002818, 0.003591, 0.003968, 0.004203, 0.001506, 0.001766, 0.016925, 0.033010, 0.058397, 0.035632, 0.011234, 0.003230, 0.009993, 0.002086, 0.006319, 0.003074, 0.003216, 0.004092, 0.004104, 0.002995, 0.002160, 0.001492, -1.584530, 0.199432, 0.444164, 0.184670, 0.002556, 0.002069, 0.003553, 0.001793, 0.005501, 0.005072, 0.003964, 0.004594, 0.001764, 0.004516, 0.000774, 0.000726, 0.016601, 0.014883, 0.029983, 0.012714, 0.005580, 0.005128, 0.005481, 0.003655, 0.001979, 0.002564, 0.002159, 0.001855, 0.002083, 0.180994, 0.020807, 0.150360, 0.016826,
                          0.001640, 0.001543, 0.004984, 0.000968, 0.002310, 0.012149, 0.003316, 0.001902, 0.000806, 0.003047, 0.001007, 0.000877, 0.008201, 0.085071, 0.049895, 0.016940, 0.002507, 0.009090, 0.011354, 0.003151, 0.002471, 0.010308, 0.001818, 0.001478, 0.000738, 0.007135, 0.001360, 0.000516, 0.059667, -1.170780, 0.335384, 0.150181, 0.001178, 0.001656, 0.003341, 0.000602, 0.002488, 0.018559, 0.003134, 0.003252, 0.000579, 0.002579, 0.000989, 0.000651, 0.004389, 0.045858, 0.027239, 0.006564, 0.006871, 0.001753, 0.000914, 0.005262, 0.001523, 0.001822, 0.004500, 0.002840, 0.001328, 0.044211, 0.059985, 0.113983, 0.010919,
                          0.000975, 0.001109, 0.005684, 0.000681, 0.001913, 0.004572, 0.003291, 0.001317, 0.000725, 0.002806, 0.002330, 0.000972, 0.009593, 0.046325, 0.065722, 0.019274, 0.003308, 0.002634, 0.018865, 0.001135, 0.001553, 0.002760, 0.002400, 0.001175, 0.001015, 0.002941, 0.003721, 0.000851, 0.074741, 0.188631, -0.904017, 0.096500, 0.001086, 0.001079, 0.002661, 0.000256, 0.002255, 0.005230, 0.007045, 0.002148, 0.000437, 0.001095, 0.000855, 0.000231, 0.003592, 0.013058, 0.054897, 0.005046, 0.004292, 0.001680, 0.001265, 0.001862, 0.002565, 0.000824, 0.002028, 0.003398, 0.000498, 0.037640, 0.021434, 0.141955, 0.014083,
                          0.007209, 0.001323, 0.005610, 0.004886, 0.003793, 0.003810, 0.002521, 0.003561, 0.002215, 0.001267, 0.001660, 0.001665, 0.008391, 0.029481, 0.041228, 0.069758, 0.010019, 0.003265, 0.008282, 0.005223, 0.001983, 0.002835, 0.001759, 0.008452, 0.001413, 0.003226, 0.001595, 0.003338, 0.088483, 0.240513, 0.274777, -1.339180, 0.003745, 0.000851, 0.002341, 0.001435, 0.005900, 0.005917, 0.002263, 0.009201, 0.000952, 0.000804, 0.000579, 0.000868, 0.006795, 0.012082, 0.015077, 0.026331, 0.003221, 0.006827, 0.001106, 0.002110, 0.000827, 0.002982, 0.000486, 0.001516, 0.001824, 0.149304, 0.016458, 0.179996, 0.033836,
                          0.052712, 0.009301, 0.023870, 0.013485, 0.008621, 0.008217, 0.003627, 0.008269, 0.004306, 0.007646, 0.001371, 0.007195, 0.002531, 0.000551, 0.003842, 0.001842, 0.051581, 0.004169, 0.023491, 0.003625, 0.007344, 0.002598, 0.003372, 0.003000, 0.001196, 0.001906, 0.001036, 0.001388, 0.000475, 0.000731, 0.001198, 0.001451, -0.909708, 0.051581, 0.414485, 0.075046, 0.012436, 0.008018, 0.005118, 0.013648, 0.011244, 0.003260, 0.002304, 0.002853, 0.004717, 0.001663, 0.004077, 0.003948, 0.001111, 0.002447, 0.008217, 0.003565, 0.003744, 0.005991, 0.000044, 0.001035, 0.000062, 0.000913, 0.000668, 0.001035, 0.000532,
                          0.008857, 0.070427, 0.018306, 0.022056, 0.002690, 0.012110, 0.004484, 0.003012, 0.000869, 0.026175, 0.001588, 0.007613, 0.000272, 0.001910, 0.001409, 0.000530, 0.006065, 0.015397, 0.012571, 0.004629, 0.003574, 0.005548, 0.003498, 0.003391, 0.000993, 0.002335, 0.000926, 0.000669, 0.000407, 0.001088, 0.001261, 0.000349, 0.054649, -0.887159, 0.113244, 0.366986, 0.003726, 0.013919, 0.003567, 0.004327, 0.006999, 0.025744, 0.007086, 0.005164, 0.000556, 0.003034, 0.001970, 0.001121, 0.002863, 0.001857, 0.002944, 0.009274, 0.003594, 0.005240, 0.000088, 0.001469, 0.000123, 0.000460, 0.001223, 0.000518, 0.000405,
                          0.012531, 0.010163, 0.043732, 0.005785, 0.003946, 0.009429, 0.005820, 0.004362, 0.001577, 0.010958, 0.003869, 0.004699, 0.000833, 0.000902, 0.006177, 0.000748, 0.017852, 0.006797, 0.076676, 0.003487, 0.004516, 0.005231, 0.004397, 0.001978, 0.001627, 0.004181, 0.004485, 0.001399, 0.000529, 0.001661, 0.002352, 0.000727, 0.332074, 0.085634, -0.827623, 0.041741, 0.005780, 0.017439, 0.013747, 0.007609, 0.002938, 0.004955, 0.006328, 0.001732, 0.002173, 0.003341, 0.010932, 0.002150, 0.001869, 0.000733, 0.003117, 0.006317, 0.005270, 0.004059, 0.000168, 0.000740, 0.000058, 0.000648, 0.000494, 0.001538, 0.000613,
                          0.012536, 0.027045, 0.010321, 0.056647, 0.003672, 0.004824, 0.003299, 0.007199, 0.002037, 0.007188, 0.000825, 0.011514, 0.000261, 0.000605, 0.001393, 0.001232, 0.012303, 0.005715, 0.009001, 0.006168, 0.005558, 0.001978, 0.001946, 0.004865, 0.000847, 0.001032, 0.000627, 0.001027, 0.000400, 0.000449, 0.000340, 0.000668, 0.090138, 0.416040, 0.062578, -0.870644, 0.003111, 0.005763, 0.002636, 0.008120, 0.011860, 0.007154, 0.004118, 0.010285, 0.000834, 0.000876, 0.001336, 0.004188, 0.001246, 0.002364, 0.007991, 0.007123, 0.003352, 0.011668, 0.000187, 0.000394, 0.000066, 0.000267, 0.000650, 0.001134, 0.001641,
                          0.017303, 0.006179, 0.013392, 0.006338, 0.022620, 0.012553, 0.007019, 0.010843, 0.003845, 0.008243, 0.002812, 0.012257, 0.004365, 0.003361, 0.012101, 0.006234, 0.014920, 0.002671, 0.009156, 0.003235, 0.015856, 0.005384, 0.006728, 0.006976, 0.002362, 0.001306, 0.001424, 0.001514, 0.002248, 0.003399, 0.005476, 0.005033, 0.027375, 0.007742, 0.015882, 0.005701, -1.281480, 0.271180, 0.121716, 0.349548, 0.025036, 0.012072, 0.009808, 0.009767, 0.019693, 0.013435, 0.033049, 0.016723, 0.001508, 0.001516, 0.033088, 0.018297, 0.011147, 0.019128, 0.002102, 0.002326, 0.004203, 0.003712, 0.002389, 0.004174, 0.002013,
                          0.004339, 0.012881, 0.011766, 0.003063, 0.004790, 0.038633, 0.009271, 0.008408, 0.001163, 0.020746, 0.001899, 0.005334, 0.001503, 0.008937, 0.013491, 0.002769, 0.003438, 0.003622, 0.012598, 0.001956, 0.002308, 0.015163, 0.004317, 0.002997, 0.001267, 0.004505, 0.002025, 0.000943, 0.000985, 0.012045, 0.006035, 0.002398, 0.008386, 0.013740, 0.022765, 0.005018, 0.128837, -1.030450, 0.130039, 0.243982, 0.006334, 0.030818, 0.005660, 0.007301, 0.003441, 0.044243, 0.025613, 0.007859, 0.005323, 0.001728, 0.009658, 0.050009, 0.013015, 0.016277, 0.008300, 0.001580, 0.001481, 0.000944, 0.008794, 0.001789, 0.001926,
                          0.009566, 0.011891, 0.031600, 0.006577, 0.014167, 0.020223, 0.023350, 0.008562, 0.001623, 0.014675, 0.005167, 0.008501, 0.002513, 0.004636, 0.021795, 0.003175, 0.006020, 0.004121, 0.029241, 0.005961, 0.008266, 0.008675, 0.019554, 0.003597, 0.001622, 0.003900, 0.003378, 0.000952, 0.002404, 0.006353, 0.025396, 0.002864, 0.016719, 0.010999, 0.056056, 0.007170, 0.180631, 0.406194, -1.458140, 0.207610, 0.008643, 0.024944, 0.009536, 0.008580, 0.008847, 0.025380, 0.041698, 0.007592, 0.003135, 0.001237, 0.015863, 0.024816, 0.022260, 0.023773, 0.006593, 0.000907, 0.002668, 0.002343, 0.004267, 0.006127, 0.003222,
                          0.011602, 0.006004, 0.014616, 0.007243, 0.009141, 0.011263, 0.004977, 0.020036, 0.002668, 0.006236, 0.001699, 0.010036, 0.001894, 0.002947, 0.013725, 0.006957, 0.008868, 0.001765, 0.007258, 0.002147, 0.005784, 0.003316, 0.002992, 0.010379, 0.000813, 0.001519, 0.001603, 0.001265, 0.001255, 0.002969, 0.003487, 0.005246, 0.020080, 0.006009, 0.013972, 0.009947, 0.233624, 0.343229, 0.093501, -1.138310, 0.008293, 0.010152, 0.005295, 0.016253, 0.006961, 0.015340, 0.017163, 0.031158, 0.002042, 0.001846, 0.019846, 0.024005, 0.009984, 0.038159, 0.003614, 0.000576, 0.004391, 0.001681, 0.002354, 0.003321, 0.003799,
                          0.007718, 0.013910, 0.010894, 0.015057, 0.003265, 0.002856, 0.001713, 0.003015, 0.006239, 0.010956, 0.002514, 0.008628, 0.000414, 0.000590, 0.004146, 0.000985, 0.007471, 0.001928, 0.006395, 0.004366, 0.002513, 0.001947, 0.001546, 0.001922, 0.001872, 0.001011, 0.001155, 0.000788, 0.000597, 0.000654, 0.000879, 0.000672, 0.020477, 0.012030, 0.006679, 0.017983, 0.020712, 0.011030, 0.004818, 0.010265, -1.001990, 0.242481, 0.193001, 0.290581, 0.000570, 0.001342, 0.001968, 0.001121, 0.001111, 0.000979, 0.007586, 0.003922, 0.003928, 0.011724, 0.001627, 0.001603, 0.001521, 0.000366, 0.001458, 0.000658, 0.001838,
                          0.003801, 0.029987, 0.012139, 0.007007, 0.001814, 0.006092, 0.001982, 0.001902, 0.001136, 0.035981, 0.001634, 0.009677, 0.000400, 0.001044, 0.004240, 0.000525, 0.002058, 0.004518, 0.007542, 0.002031, 0.001823, 0.004818, 0.003360, 0.001966, 0.001416, 0.003403, 0.001947, 0.001064, 0.001138, 0.002173, 0.001640, 0.000423, 0.004426, 0.032982, 0.008394, 0.008085, 0.007444, 0.039997, 0.010364, 0.009366, 0.180725, -0.893499, 0.162152, 0.220989, 0.000952, 0.004233, 0.004363, 0.001619, 0.001751, 0.001882, 0.002320, 0.010545, 0.004139, 0.005933, 0.001842, 0.002790, 0.000292, 0.000318, 0.003491, 0.000432, 0.000992,
                          0.002099, 0.021315, 0.011383, 0.012858, 0.002607, 0.005139, 0.003528, 0.003445, 0.004156, 0.023086, 0.012829, 0.012873, 0.000424, 0.000968, 0.006831, 0.001134, 0.004042, 0.004606, 0.016430, 0.005227, 0.002779, 0.003677, 0.004098, 0.003590, 0.001782, 0.002418, 0.005217, 0.001421, 0.000421, 0.001800, 0.002767, 0.000658, 0.006759, 0.019620, 0.023171, 0.010059, 0.013071, 0.015877, 0.008564, 0.010558, 0.310904, 0.350466, -1.194700, 0.187362, 0.001220, 0.002572, 0.012307, 0.001980, 0.001040, 0.001711, 0.003135, 0.006574, 0.004036, 0.006286, 0.000900, 0.005809, 0.000218, 0.000354, 0.001905, 0.001128, 0.001507,
                          0.003913, 0.010936, 0.004503, 0.017684, 0.002025, 0.001950, 0.000642, 0.003265, 0.001330, 0.010388, 0.001203, 0.011430, 0.000132, 0.000493, 0.001255, 0.000630, 0.003793, 0.001103, 0.001564, 0.001311, 0.001467, 0.000998, 0.000474, 0.001358, 0.001005, 0.000951, 0.000616, 0.000943, 0.000248, 0.000744, 0.000470, 0.000620, 0.005256, 0.008979, 0.003983, 0.015775, 0.008173, 0.012861, 0.004838, 0.020351, 0.293935, 0.299925, 0.117652, -0.910018, 0.000359, 0.002136, 0.001832, 0.001882, 0.000719, 0.001709, 0.002387, 0.002737, 0.002236, 0.008580, 0.000484, 0.000704, 0.001019, 0.000314, 0.000642, 0.000455, 0.000651,
                          0.010467, 0.004056, 0.010637, 0.002905, 0.021502, 0.017523, 0.007659, 0.013462, 0.004259, 0.003591, 0.001820, 0.002015, 0.073626, 0.045786, 0.025733, 0.055936, 0.007556, 0.003008, 0.005937, 0.003408, 0.002580, 0.003293, 0.002172, 0.001234, 0.002211, 0.001796, 0.001173, 0.001584, 0.011640, 0.010286, 0.014968, 0.009944, 0.017814, 0.001984, 0.010243, 0.002624, 0.033787, 0.012424, 0.010227, 0.017869, 0.001183, 0.002650, 0.001570, 0.000735, -1.493460, 0.203898, 0.345114, 0.368061, 0.002358, 0.003131, 0.004999, 0.004287, 0.003022, 0.003004, 0.003421, 0.002632, 0.011550, 0.013059, 0.008267, 0.014380, 0.005397,
                          0.005405, 0.002791, 0.008583, 0.001164, 0.007552, 0.039220, 0.008558, 0.006143, 0.000940, 0.005648, 0.000828, 0.001193, 0.015323, 0.174351, 0.017427, 0.044926, 0.002741, 0.004883, 0.003379, 0.001549, 0.000876, 0.004402, 0.001720, 0.001467, 0.000632, 0.004311, 0.000589, 0.000815, 0.004033, 0.041533, 0.021028, 0.006833, 0.002427, 0.004179, 0.006087, 0.001065, 0.008907, 0.061739, 0.011338, 0.015217, 0.001076, 0.004552, 0.001280, 0.001692, 0.078794, -1.212990, 0.303670, 0.203144, 0.007099, 0.002857, 0.001471, 0.004969, 0.001776, 0.001700, 0.006424, 0.001449, 0.003064, 0.002814, 0.019319, 0.009426, 0.004614,
                          0.002595, 0.001966, 0.007534, 0.001581, 0.008624, 0.014124, 0.013324, 0.005337, 0.001035, 0.003351, 0.003186, 0.001286, 0.023173, 0.085029, 0.046076, 0.042407, 0.003561, 0.002891, 0.007949, 0.002581, 0.000959, 0.002048, 0.002193, 0.001242, 0.000710, 0.001513, 0.002685, 0.000730, 0.005957, 0.018089, 0.064818, 0.006252, 0.004363, 0.001990, 0.014601, 0.001190, 0.016066, 0.026207, 0.013659, 0.012483, 0.001156, 0.003439, 0.004489, 0.001064, 0.097787, 0.222658, -0.992739, 0.131449, 0.004976, 0.002213, 0.000945, 0.003502, 0.002068, 0.002101, 0.002046, 0.001658, 0.000851, 0.003054, 0.011043, 0.014682, 0.006192,
                          0.005572, 0.002577, 0.005992, 0.004873, 0.009758, 0.011682, 0.005699, 0.020885, 0.001672, 0.002882, 0.000651, 0.001340, 0.026730, 0.060459, 0.016867, 0.138357, 0.004548, 0.001129, 0.002439, 0.001458, 0.002612, 0.002247, 0.000990, 0.002160, 0.000599, 0.001149, 0.000616, 0.000940, 0.004283, 0.007391, 0.010101, 0.018512, 0.007163, 0.001919, 0.004868, 0.006327, 0.013783, 0.013634, 0.004217, 0.038423, 0.001116, 0.002164, 0.001224, 0.001854, 0.176821, 0.252543, 0.222870, -1.191070, 0.001942, 0.002855, 0.002745, 0.002012, 0.002085, 0.004058, 0.003119, 0.001003, 0.004932, 0.008134, 0.007219, 0.010875, 0.013991,
                          0.001200, 0.009012, 0.003402, 0.003111, 0.001477, 0.006545, 0.002898, 0.001601, 0.000823, 0.005249, 0.001227, 0.000855, 0.001095, 0.005034, 0.008702, 0.001471, 0.002123, 0.040938, 0.006174, 0.011699, 0.000853, 0.001218, 0.000345, 0.000926, 0.000613, 0.005607, 0.000897, 0.001605, 0.001350, 0.005556, 0.006171, 0.001627, 0.001448, 0.003521, 0.003039, 0.001352, 0.000893, 0.006632, 0.001250, 0.001809, 0.000795, 0.001680, 0.000462, 0.000509, 0.000814, 0.006338, 0.006059, 0.001394, -0.644620, 0.302981, 0.001357, 0.007620, 0.001670, 0.002468, 0.004368, 0.015225, 0.000389, 0.000956, 0.100380, 0.002562, 0.025244,
                          0.004065, 0.003457, 0.006518, 0.008629, 0.003710, 0.003118, 0.002812, 0.002363, 0.001588, 0.002712, 0.001987, 0.004431, 0.002587, 0.004053, 0.005988, 0.004384, 0.006010, 0.021832, 0.005496, 0.036822, 0.000868, 0.001304, 0.000895, 0.002285, 0.000515, 0.001632, 0.001142, 0.001780, 0.002056, 0.002350, 0.004004, 0.005713, 0.005285, 0.003786, 0.001976, 0.004250, 0.001487, 0.003567, 0.000818, 0.002709, 0.001161, 0.002995, 0.001259, 0.002003, 0.001790, 0.004227, 0.004466, 0.003398, 0.502104, -0.889242, 0.003720, 0.004493, 0.002479, 0.009767, 0.001267, 0.018544, 0.004176, 0.003871, 0.047826, 0.002770, 0.085944,
                          0.022964, 0.014775, 0.013255, 0.018534, 0.072694, 0.021831, 0.014502, 0.036591, 0.005888, 0.041008, 0.001955, 0.049888, 0.002873, 0.001379, 0.009188, 0.001824, 0.016817, 0.004762, 0.012932, 0.007409, 0.033163, 0.006954, 0.006778, 0.007195, 0.004420, 0.002458, 0.001675, 0.002394, 0.003144, 0.001753, 0.004313, 0.001324, 0.025381, 0.008583, 0.012017, 0.020550, 0.046431, 0.028525, 0.014999, 0.041667, 0.012867, 0.005280, 0.003301, 0.004002, 0.004089, 0.003114, 0.002728, 0.004673, 0.003216, 0.005322, -1.246680, 0.192909, 0.098011, 0.239222, 0.001093, 0.001864, 0.002915, 0.005421, 0.003352, 0.006357, 0.002154,
                          0.006451, 0.021710, 0.015353, 0.007624, 0.019121, 0.078194, 0.015972, 0.018985, 0.001075, 0.084504, 0.002187, 0.029025, 0.001338, 0.002990, 0.009145, 0.001011, 0.010410, 0.011286, 0.016494, 0.004324, 0.005078, 0.024276, 0.005134, 0.004161, 0.001647, 0.008802, 0.001895, 0.001039, 0.001243, 0.005980, 0.003763, 0.001498, 0.006529, 0.016030, 0.014439, 0.010861, 0.015221, 0.087565, 0.013911, 0.029879, 0.003944, 0.014226, 0.004104, 0.002721, 0.002079, 0.006236, 0.005994, 0.002031, 0.010709, 0.003810, 0.114363, -1.153740, 0.131575, 0.202385, 0.006684, 0.004617, 0.002342, 0.000649, 0.014746, 0.001673, 0.002707,
                          0.007082, 0.015432, 0.032219, 0.012789, 0.033765, 0.055041, 0.048284, 0.021615, 0.002700, 0.094500, 0.004110, 0.056467, 0.001681, 0.002281, 0.011038, 0.001689, 0.010734, 0.010575, 0.030719, 0.004774, 0.009082, 0.011108, 0.014826, 0.009788, 0.003299, 0.008087, 0.006113, 0.003097, 0.001480, 0.003806, 0.011401, 0.001291, 0.015080, 0.013661, 0.026493, 0.011241, 0.020395, 0.050122, 0.027444, 0.027333, 0.008686, 0.012283, 0.005541, 0.004889, 0.003222, 0.004900, 0.007782, 0.004629, 0.005161, 0.004624, 0.127796, 0.289390, -1.480360, 0.239634, 0.005493, 0.005746, 0.001772, 0.001551, 0.005466, 0.004601, 0.004546,
                          0.012665, 0.012137, 0.015248, 0.018067, 0.024096, 0.029488, 0.013503, 0.048101, 0.002453, 0.046646, 0.001570, 0.055765, 0.001583, 0.001851, 0.004513, 0.003210, 0.016172, 0.006506, 0.009553, 0.009362, 0.010711, 0.006902, 0.004384, 0.017472, 0.002480, 0.002938, 0.000999, 0.003730, 0.001124, 0.002670, 0.002147, 0.002729, 0.014149, 0.011680, 0.011963, 0.022941, 0.020521, 0.036754, 0.017185, 0.061251, 0.015203, 0.010322, 0.005060, 0.011000, 0.001879, 0.002751, 0.004636, 0.005281, 0.004473, 0.010682, 0.182890, 0.260996, 0.140506, -1.271040, 0.001663, 0.001754, 0.004489, 0.001936, 0.006448, 0.001867, 0.009983,
                          0.000176, 0.001824, 0.000421, 0.000858, 0.000698, 0.002653, 0.001888, 0.000667, 0.000397, 0.010403, 0.000667, 0.001471, 0.000608, 0.004724, 0.005471, 0.001222, 0.000244, 0.000998, 0.000407, 0.000846, 0.000203, 0.000536, 0.000437, 0.000097, 0.000171, 0.001672, 0.001350, 0.000211, 0.000585, 0.004073, 0.003265, 0.000274, 0.000064, 0.000121, 0.000306, 0.000227, 0.001393, 0.011576, 0.002944, 0.003583, 0.001303, 0.001980, 0.000448, 0.000383, 0.001321, 0.006420, 0.002789, 0.002507, 0.004890, 0.000856, 0.000516, 0.005324, 0.001989, 0.001027, -0.466903, 0.000873, 0.357340, 0.000823, 0.005093, 0.000286, 0.001006,
                          0.002247, 0.001360, 0.001324, 0.000608, 0.000469, 0.002072, 0.001264, 0.001555, 0.000784, 0.001717, 0.005216, 0.000676, 0.000762, 0.002136, 0.009796, 0.000923, 0.003962, 0.004638, 0.005973, 0.003629, 0.000839, 0.000486, 0.001086, 0.000908, 0.000632, 0.003437, 0.004729, 0.000209, 0.000629, 0.003220, 0.006850, 0.001073, 0.001891, 0.002532, 0.001688, 0.000599, 0.001931, 0.002759, 0.000507, 0.000715, 0.001608, 0.003755, 0.003617, 0.000698, 0.001273, 0.001814, 0.002831, 0.001010, 0.021345, 0.015687, 0.001102, 0.004606, 0.002606, 0.001357, 0.001093, -0.196689, 0.001299, 0.002416, 0.019373, 0.007356, 0.014010,
                          0.000190, 0.001600, 0.000190, 0.001653, 0.001322, 0.000941, 0.000872, 0.002271, 0.000333, 0.002470, 0.000195, 0.003956, 0.000456, 0.002396, 0.002018, 0.001832, 0.000808, 0.000914, 0.000817, 0.000386, 0.000153, 0.000111, 0.000239, 0.000224, 0.000267, 0.000182, 0.000610, 0.001051, 0.000777, 0.001654, 0.001104, 0.001419, 0.000124, 0.000234, 0.000144, 0.000110, 0.003833, 0.002843, 0.001640, 0.005992, 0.001677, 0.000433, 0.000149, 0.001110, 0.006141, 0.004216, 0.001596, 0.005458, 0.000599, 0.003883, 0.001895, 0.002568, 0.000883, 0.003817, 0.491875, 0.001428, -0.582349, 0.001114, 0.002527, 0.001361, 0.001287,
                          0.009698, 0.002434, 0.004231, 0.002181, 0.007355, 0.003706, 0.001486, 0.003563, 0.003005, 0.001332, 0.000938, 0.000850, 0.036737, 0.022514, 0.050124, 0.039565, 0.011123, 0.001684, 0.005324, 0.003985, 0.002657, 0.000688, 0.001101, 0.001412, 0.000924, 0.001225, 0.001344, 0.001364, 0.134940, 0.110170, 0.166768, 0.232319, 0.003665, 0.001742, 0.003248, 0.000893, 0.006771, 0.003626, 0.002881, 0.004588, 0.000807, 0.000941, 0.000485, 0.000684, 0.013885, 0.007741, 0.011459, 0.018003, 0.002946, 0.007198, 0.007047, 0.001422, 0.001547, 0.003293, 0.002265, 0.005311, 0.002229, -1.350170, 0.022013, 0.315940, 0.030797,
                          0.000913, 0.002054, 0.001345, 0.000711, 0.002252, 0.006692, 0.002769, 0.001385, 0.000329, 0.002140, 0.000747, 0.000745, 0.003139, 0.026339, 0.017464, 0.005509, 0.000670, 0.010085, 0.002403, 0.002351, 0.000292, 0.002068, 0.000435, 0.000435, 0.000351, 0.001420, 0.000803, 0.000370, 0.004178, 0.040259, 0.025578, 0.006897, 0.000722, 0.001248, 0.000668, 0.000585, 0.001174, 0.009094, 0.001413, 0.001731, 0.000866, 0.002782, 0.000702, 0.000377, 0.002367, 0.014317, 0.011161, 0.004303, 0.083316, 0.023953, 0.001174, 0.008709, 0.001468, 0.002953, 0.003776, 0.011470, 0.001361, 0.005929, -0.674727, 0.014580, 0.289399,
                          0.005134, 0.001612, 0.009549, 0.003780, 0.007256, 0.004386, 0.003892, 0.003650, 0.002703, 0.002099, 0.002294, 0.001057, 0.014394, 0.026127, 0.070423, 0.029953, 0.004552, 0.001303, 0.007416, 0.001937, 0.001483, 0.001460, 0.001658, 0.000831, 0.001244, 0.002063, 0.001943, 0.001452, 0.062751, 0.158996, 0.352069, 0.156779, 0.002326, 0.001100, 0.004314, 0.002122, 0.004262, 0.003845, 0.004216, 0.005074, 0.000813, 0.000716, 0.000864, 0.000555, 0.008558, 0.014517, 0.030840, 0.013473, 0.004420, 0.002884, 0.004626, 0.002053, 0.002568, 0.001777, 0.000441, 0.009051, 0.001523, 0.176854, 0.030303, -1.304410, 0.024070,
                          0.002251, 0.002443, 0.002279, 0.003053, 0.002903, 0.003136, 0.001739, 0.002095, 0.000480, 0.002546, 0.000343, 0.001933, 0.004630, 0.012038, 0.017844, 0.022139, 0.001431, 0.004439, 0.001592, 0.006699, 0.000860, 0.000519, 0.001155, 0.001232, 0.000324, 0.000767, 0.000553, 0.000376, 0.005670, 0.012297, 0.028200, 0.023794, 0.000966, 0.000694, 0.001388, 0.002479, 0.001660, 0.003343, 0.001790, 0.004686, 0.001831, 0.001327, 0.000932, 0.000642, 0.002593, 0.005737, 0.010501, 0.013995, 0.035159, 0.072228, 0.001266, 0.002683, 0.002048, 0.007672, 0.001252, 0.013918, 0.001164, 0.013919, 0.485610, 0.019433, -0.882676 };

    Db_matrix *charQ  = new Db_matrix(char_as,char_as,"Q_char");

    for(int j=0;j<char_as;j++)
    {
        for(int i=0;i<char_as;i++)
        {
            charQ->s(tmp_q[j*char_as+i],j,i);
        }
    }

    // Find eigenvalues and eigenvectors.
    //
    charU = new Db_matrix(char_as,char_as,"eigenvectors_1");
    charU->initialise();
    charV = new Db_matrix(char_as,char_as,"eigenvectors_2");
    charV->initialise();
    charRoot = new Db_matrix(char_as,"eigenvalues");
    charRoot->initialise();

    build_model(char_as,charPi,charQ,charU,charV,charRoot);

    if (Settings::noise > 4) {
        print_char_q_matrices(charQ);
    }

    delete charQ;

    character_alphabet = this->get_codon_character_alphabet();
    full_character_alphabet = this->get_codon_full_character_alphabet();

    ancestral_character_alphabet.clear();

    std::string full_alpha = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTTNNN";

    for(int i=0;i<62;i++)
    {
        ancestral_character_alphabet.push_back(full_alpha.substr(i*3,3));
    }

//    for(int i=0;i<61;i++)
//        for(int j=0;j<61;j++)
//            if(charPi->g(i) > charPi->g(j))
//                ancestral_character_alphabet.push_back(full_alpha.substr(i*3,3));
//            else
//                ancestral_character_alphabet.push_back(full_alpha.substr(j*3,3));


    string alpha = "xACxGxxxT";
    string ambiguity = "xxAMCRSxGWYxKxxxT";

    for(int i=0;i<60;i++)
    {
        for(int j=i+1;j<61;j++)
        {
            string codon;

            string c1 = full_alpha.substr(i*3,1);
            string c2 = full_alpha.substr(j*3,1);

            string::size_type loc1 = alpha.find(c1,0);
            string::size_type loc2 = alpha.find(c2,0);

            if(loc1 != string::npos && loc2 != string::npos)
                codon += ambiguity.at(int(loc1+loc2));
            else
                codon += "N";


            c1 = full_alpha.substr(i*3+1,1);
            c2 = full_alpha.substr(j*3+1,1);

            loc1 = alpha.find(c1,0);
            loc2 = alpha.find(c2,0);

            if(loc1 != string::npos && loc2 != string::npos)
                codon += ambiguity.at(int(loc1+loc2));
            else
                codon += "N";


            c1 = full_alpha.substr(i*3+2,1);
            c2 = full_alpha.substr(j*3+2,1);

            loc1 = alpha.find(c1,0);
            loc2 = alpha.find(c2,0);

            if(loc1 != string::npos && loc2 != string::npos)
                codon += ambiguity.at(int(loc1+loc2));
            else
                codon += "N";

            ancestral_character_alphabet.push_back(codon);

        }
    }

}

/*******************************************/

void Model_factory::build_model(int s,Db_matrix *pi,Db_matrix *q,Db_matrix *wU,Db_matrix *wV,Db_matrix *wRoot)
{


    // Find eigenvalues and eigenvectors.
    //

    Eigen* e = new Eigen();

    double* tpi = new double[s];
    double* tsq = new double[s];
    double* tcq = new double[s*s];
    double* twr = new double[s];
    double* twu = new double[s*s];
    double* twv = new double[s*s];
    int npi0 = 0;

    for(int i=0;i<s;i++)
    {
        tpi[i] = pi->g(i);

        for(int j=0;j<s;j++)
        {
            tcq[i*s+j] = q->g(i,j);
        }
    }

    if (e->getpi_sqrt (tpi, tsq, s, &npi0) != 0) {
        Log_output::write_out("Model_factory: Error in eigen square roots!!\n",0);
        exit (1);
    }


    if (e->eigenQREV (tcq, tpi, tsq, s, npi0, twr, twu, twv) != 0) {
        Log_output::write_out("Model_factory: Error in eigen QREV!!\n",0);
        exit (1);
    }

    for(int i=0;i<s;i++)
    {
        wRoot->s(twr[i],i);

        for(int j=0;j<s;j++)
        {
            wU->s(twu[i*s+j],i,j);
            wV->s(twv[i*s+j],i,j);
        }
    }

    delete []tpi;
    delete []tsq;
    delete []tcq;
    delete []twr;
    delete []twu;
    delete []twv;

    delete e;

}

/*******************************************/

Evol_model Model_factory::alignment_model(double distance)
{

    // Compute the P matrix for regular DNA alphabet (four bases).
    //
    Eigen* e = new Eigen();

    double* tmr = new double[char_as*char_as];
    double* twr = new double[char_as];
    double* twu = new double[char_as*char_as];
    double* twv = new double[char_as*char_as];

    for(int i=0;i<char_as;i++)
    {
        twr[i] = charRoot->g(i);

        for(int j=0;j<char_as;j++)
        {
            twu[i*char_as+j] = charU->g(i,j);
            twv[i*char_as+j] = charV->g(i,j);
        }
    }

    e->computePMatrix(char_as,tmr,twu,twv,twr,distance);

    Evol_model model(sequence_data_type, distance);

    model.log_ext_prob = log(char_ext_prob);
    model.ext_prob = char_ext_prob;

    if( (Settings_handle::st.is("454") || Settings_handle::st.is("homopolymer"))&& Settings_handle::st.is("pileup-alignment"))
    {
        char_ins_rate = 0.25;
        char_del_rate = 0.25;
    }

    model.ins_rate = char_ins_rate;
    model.del_rate = char_del_rate;

    model.ins_prob = (1.0-exp(-1.0*char_ins_rate*distance));
    model.del_prob = (1.0-exp(-1.0*char_del_rate*distance));

    double t = (1.0-exp(-0.5*(char_ins_rate+char_del_rate)*distance));

    model.log_id_prob = log(t);
    model.log_match_prob = log(1.0-2*t);
    model.id_prob = t;
    model.match_prob = 1.0-2*t;


    model.log_end_ext_prob = log(char_end_ext_prob);
    model.log_break_ext_prob = log(char_break_ext_prob);

    model.end_ext_prob = char_end_ext_prob;
    model.break_ext_prob = char_break_ext_prob;


    for(int i=0;i<char_as;i++)
    {
        model.charPi->s(charPi->g(i),i);
        model.logCharPi->s(log( charPi->g(i) ),i);

        for(int j=0;j<char_as;j++)
        {
            model.mostcommon_table->s(mostcommon_table->g(i,j),i,j);

            float sp = tmr[i*char_as+j];
            if( Settings_handle::st.is("no-score-scaling") )
            {
                float lo = sp / ( charPi->g(i) * charPi->g(j) );
                model.charPr->s(lo,i,j);
                model.logCharPr->s(log( lo ),i,j);
            }
            else if( ! Settings_handle::st.is("no-log-odds") )
            {
                float lo = 0.5 * ( charPi->g(i) + charPi->g(j) ) * sp / ( charPi->g(i) * charPi->g(j) );
                model.charPr->s(lo,i,j);
                model.logCharPr->s(log( lo ),i,j);
            }
            else
            {
                model.charPr->s(tmr[i*char_as+j],i,j);
                model.logCharPr->s(log( tmr[i*char_as+j] ),i,j);
            }
        }
    }

    delete[] tmr;
    delete[] twr;
    delete[] twu;
    delete[] twv;

    delete e;

    /***************************************************************/

    if(sequence_data_type == Model_factory::dna && !Settings_handle::st.is("codons"))
    {
        // Extend the P matrix for ambiguity DNA alphabet (16 bases).
        //

        float ambiguity_factor = 1.0;
        if( Settings_handle::st.is("ambiguity-factor") )
            ambiguity_factor = Settings_handle::st.get("ambiguity-factor").as<float>();

        if( ambiguity_factor > 1.0 || ambiguity_factor < 0 )
            ambiguity_factor = 1.0;

        for(unsigned int ai=0;ai<char_symbols.size();ai++)
        {
            Char_symbol *a = &char_symbols.at(ai);

            float probability = pow(ambiguity_factor, a->n_units);

            for(int aj=0;aj<a->n_units;aj++){
                int at = char_alphabet.find(a->residues.at(aj));
                char_ambiguity->s(probability, at, ai);
            }
        }



        for(int i=0;i<char_fas;i++)
        {
            for(int j=0;j<char_fas;j++)
            {
                model.parsimony_table->s( parsimony_table->g(i,j), i,j);

                if(i<char_as && j<char_as)
                    continue;

                double max = 0;

                for(int n=0;n<char_as;n++)
                {
                    for(int m=0;m<char_as;m++)
                    {
                        double t = model.charPr->g(n,m)*char_ambiguity->g(m,j)*char_ambiguity->g(n,i);
                        if(max < t) max = t;
                    }
                }

                model.charPr->s(max,i,j);
                model.logCharPr->s(log(max),i,j);
            }
        }

    }

    /***************************************************************/

    else if( ( sequence_data_type == Model_factory::dna && Settings_handle::st.is("codons") ) ||
            sequence_data_type == Model_factory::codon )
    {

        for(int i=0;i<char_fas;i++)
        {
            for(int j=0;j<char_fas;j++)
            {
                model.parsimony_table->s( parsimony_table->g(i,j), i,j);

                if(i<char_as && j<char_as)
                    continue;

                double max = 0;
                double t;

                if(i==char_as) // NNN
                {
                    for(int n=0;n<char_as;n++)
                    {
                        double t = model.charPr->g(n,j);
                        if(max < t) max = t;
                    }
                }
                else if(j==char_as) // NNN
                {
                    for(int m=0;m<char_as;m++)
                    {
                        double t = model.charPr->g(i,m);
                        if(max < t) max = t;
                    }
                }
                else
                {
                    Codon_symbol *l1 = &codon_symbols.at(i);
                    Codon_symbol *l2 = &codon_symbols.at(j);

                    if(l1->n_units==1 && l2->n_units==2)
                    {
                        max = model.charPr->g(l1->first_codon,l2->first_codon);
                        t = model.charPr->g(l1->first_codon,l2->second_codon);
                        if(max < t) max = t;
                    }
                    else if(l1->n_units==2 && l2->n_units==1)
                    {
                        max = model.charPr->g(l1->first_codon,l2->first_codon);
                        t = model.charPr->g(l1->second_codon,l2->first_codon);
                        if(max < t) max = t;
                    }
                    else if(l1->n_units==2 && l2->n_units==2)
                    {
                        max = model.charPr->g(l1->first_codon,l2->first_codon);
                        t = model.charPr->g(l1->first_codon,l2->second_codon);
                        if(max < t) max = t;
                        t = model.charPr->g(l1->second_codon,l2->first_codon);
                        if(max < t) max = t;
                        t = model.charPr->g(l1->second_codon,l2->second_codon);
                        if(max < t) max = t;
                    }
                    else
                    {
                        Log_output::write_out("Model_factory: errors! "+Log_output::itos(i)+" "+Log_output::itos(j)+"\n",2);
                    }
                }

                model.charPr->s(max,i,j);
                model.logCharPr->s(log(max),i,j);
            }
        }

    }

    /***************************************************************/

    else if(sequence_data_type == Model_factory::protein)
    {
        // Large aa groups; easiest to do using ambiguity code

        if(Settings_handle::st.is("use-aa-groups"))
        {
            // Extend the P matrix for ambiguity protein alphabet (currently 51 residues).
            //

            float ambiguity_factor = 1.0;
            if( Settings_handle::st.is("ambiguity-factor") )
                ambiguity_factor = Settings_handle::st.get("ambiguity-factor").as<float>();

            if( ambiguity_factor > 1.0 || ambiguity_factor < 0 )
                ambiguity_factor = 1.0;

            for(unsigned int ai=0;ai<char_symbols.size();ai++)
            {
                Char_symbol *a = &char_symbols.at(ai);

                float probability = pow(ambiguity_factor, a->n_units);

                for(int aj=0;aj<a->n_units;aj++){
                    int at = char_alphabet.find(a->residues.at(aj));
                    char_ambiguity->s(probability, at, ai);
                }
            }


            for(int i=0;i<char_fas;i++)
            {
                for(int j=0;j<char_fas;j++)
                {
                    model.parsimony_table->s( parsimony_table->g(i,j), i,j);

                    if(i<char_as && j<char_as)
                        continue;

                    double max = 0;

                    for(int n=0;n<char_as;n++)
                    {
                        for(int m=0;m<char_as;m++)
                        {
                            double t = model.charPr->g(n,m)*char_ambiguity->g(m,j)*char_ambiguity->g(n,i);
                            if(max < t) max = t;
                        }
                    }

                    model.charPr->s(max,i,j);
                    model.logCharPr->s(log(max),i,j);
                }
            }
        }

        // Large alphabet; faster to do by exploiting the small group size

        else
        {
            for(int i=0;i<char_fas;i++)
            {
                for(int j=0;j<char_fas;j++)
                {
                    model.parsimony_table->s( parsimony_table->g(i,j), i,j);

                    if(i<char_as && j<char_as)
                        continue;

                    double max = 0;
                    double t;

                    if(i==char_as) // X
                    {
                        for(int n=0;n<char_as;n++)
                        {
                            double t = model.charPr->g(n,j);
                            if(max < t) max = t;
                        }
                    }
                    else if(j==char_as) //X
                    {
                        for(int m=0;m<char_as;m++)
                        {
                            double t = model.charPr->g(i,m);
                            if(max < t) max = t;
                        }
                    }
                    else
                    {
                        Char_symbol *l1 = &char_symbols.at(i);
                        Char_symbol *l2 = &char_symbols.at(j);

                        if(l1->n_units==1 && l2->n_units==2)
                        {
                            max = model.charPr->g(l1->first_residue,l2->first_residue);
                            t = model.charPr->g(l1->first_residue,l2->second_residue);
                            if(max < t) max = t;
                        }
                        else if(l1->n_units==2 && l2->n_units==1)
                        {
                            max = model.charPr->g(l1->first_residue,l2->first_residue);
                            t = model.charPr->g(l1->second_residue,l2->first_residue);
                            if(max < t) max = t;
                        }
                        else if(l1->n_units==2 && l2->n_units==2)
                        {
                            max = model.charPr->g(l1->first_residue,l2->first_residue);
                            t = model.charPr->g(l1->first_residue,l2->second_residue);
                            if(max < t) max = t;
                            t = model.charPr->g(l1->second_residue,l2->first_residue);
                            if(max < t) max = t;
                            t = model.charPr->g(l1->second_residue,l2->second_residue);
                            if(max < t) max = t;
                        }
                        else
                        {
                            Log_output::write_out("Model_factory: errors! "+Log_output::itos(i)+" "+Log_output::itos(j)+"\n",2);
                        }
                    }

                    model.charPr->s(max,i,j);
                    model.logCharPr->s(log(max),i,j);
                }
            }
        }
    }

    /***************************************************************/

    if (Settings::noise > 4) {
        print_char_p_matrices(model);
    }

    return model;
}

/*******************************************/

void Model_factory::print_int_matrix(Int_matrix *m)
{
    stringstream ss;
    ss<<"    ";
    for(int i=0;i<m->X();i++)
        ss<<" "<<setw(4)<<i;
    ss<<endl;
    for(int j=0;j<m->X();j++)
    {
        ss<<setw(4)<<j;
        for(int i=0;i<m->X();i++)
            ss<<" "<<setw(4)<<m->g(i,j);

        ss<<endl;
    }
    ss<<endl;

    Log_output::write_out(ss.str(),0);
}

void Model_factory::print_char_p_matrices(Evol_model &model)
{
    stringstream ss;
    // Print out the model
    ss<<"\nModel_factory::print_char_p_matrices()\n\n";
    ss<<"alphabet "<<char_alphabet<<endl;
    ss<<"distance "<<model.distance<<endl;

    ss<<"\ncharacter equilibrium frequencies (pi)"<<endl;
    ss << fixed << noshowpos << setprecision (4);

    for(int i=0;i<char_as;i++)
        ss << setw(4) <<char_alphabet.at(i)<<"   ";
    ss<<endl;

    for(int j=0;j<char_as;j++)
        ss<<" "<<model.charPi->g(j);
    ss<<endl;

    ss<<"\nsubstitution matrix"<<endl;
    ss << fixed << noshowpos << setprecision (4);

    for(int i=0;i<char_as;i++)
        ss << setw(6) <<char_alphabet.at(i)<<" ";
    ss<<endl;

    for(int i=0;i<char_as;i++)
    {
        ss<<char_alphabet.at(i)<<" ";
        for(int j=0;j<char_as;j++)
        {
            ss<<" "<<model.charPr->g(i,j);
        }
        ss<<endl;
    }
    ss<<"\nlog substitution matrix"<<endl;
    ss << fixed << showpos << setprecision (4);

    for(int i=0;i<char_as;i++)
        ss << setw(7) <<char_alphabet.at(i)<<" ";
    ss<<endl;

    for(int i=0;i<char_as;i++)
    {
        ss<<char_alphabet.at(i)<<" ";
        for(int j=0;j<char_as;j++)
        {
            ss<<" "<<model.logCharPr->g(i,j);
        }
        ss<<endl;
    }


    ss << noshowpos <<"\nfull alphabet "<<full_char_alphabet<<endl;

    ss<<"\nsubstitution matrix"<<endl;
    ss << fixed << noshowpos << setprecision (4);

    for(int i=0;i<char_fas;i++)
        ss << setw(6) <<full_char_alphabet.at(i)<<" ";
    ss<<endl;

    for(int i=0;i<char_fas;i++)
    {
        ss<<full_char_alphabet.at(i)<<" ";
        for(int j=0;j<char_fas;j++)
        {
            ss<<" "<<model.charPr->g(i,j);
        }
        ss<<endl;
    }


    ss<<"\nlog substitution matrix"<<endl;
    ss << fixed << showpos << setprecision (4);

    for(int i=0;i<char_fas;i++) 
        ss << setw(7) <<full_char_alphabet.at(i)<<" ";
    ss<<endl;

    for(int i=0;i<char_fas;i++)
    {
        ss<<full_char_alphabet.at(i)<<" ";
        for(int j=0;j<char_fas;j++)
        {
            ss<<" "<<model.logCharPr->g(i,j);
        }
        ss<<endl;
    }

    ss<<endl;
    ss<<"indel prob:     "<<model.id_prob<<", "<<model.log_id_prob<<endl;
    ss<<"extension prob: "<<model.ext_prob<<", "<<model.log_ext_prob<<endl;
    ss<<"end extension prob: "<<model.end_ext_prob<<", "<<model.log_end_ext_prob<<endl;
    ss<<"match prob:     "<<model.match_prob<<", "<<model.log_match_prob<<endl;;
    ss<<endl;

    Log_output::write_out(ss.str(),0);

}

void Model_factory::print_char_q_matrices(Db_matrix *charQ)
{
    stringstream ss;

    // Print out the model
    ss<<"\nModel_factory::print_char_q_matrices()\n\n";
    ss<<"alphabet "<<char_alphabet<<endl;

    ss<<"\ncharacter equilibrium frequencies (pi)"<<endl;
    ss << fixed << noshowpos << setprecision (4);

    for(int i=0;i<char_as;i++)
        ss << setw(4) <<char_alphabet.at(i)<<"   ";
    ss<<endl;

    for(int j=0;j<char_as;j++)
        ss<<" "<<charPi->g(j);
    ss<<endl;

    ss<<"\noriginal substitution matrix"<<endl;
    ss << fixed << showpos << setprecision (4);

    for(int i=0;i<char_as;i++) 
        ss << setw(7) <<char_alphabet.at(i)<<" ";
    ss<<endl;

    for(int i=0;i<char_as;i++)
    {
        ss<<full_char_alphabet.at(i)<<" ";
        for(int j=0;j<char_as;j++)
            ss<<" "<<charQ->g(i,j);

        ss<<endl;
    }

    ss<<"\neigen values & vectors"<<endl;
    ss << fixed << showpos << setprecision (4);

    for(int j=0;j<char_as;j++)
        ss<<" "<<charRoot->g(j);

    ss<<endl<<endl;

    for(int i=0;i<char_as;i++)
    {
        for(int j=0;j<char_as;j++)
            ss<<" "<<charU->g(i,j);

        ss<<endl;
    }
    ss<<endl;

    for(int i=0;i<char_as;i++)
    {
        for(int j=0;j<char_as;j++)
            ss<<" "<<charV->g(i,j);

        ss<<endl;
    }
    ss<<endl;

    ss<<"ins rate:     "<<char_ins_rate<<endl;
    ss<<"del rate:     "<<char_del_rate<<endl;
    ss<<"extension prob: "<<char_ext_prob<<endl;
    ss<<endl;

    Log_output::write_out(ss.str(),0);
}



