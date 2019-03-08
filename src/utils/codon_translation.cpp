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

#include "utils/codon_translation.h"

using namespace std;
using namespace ppa;

Codon_translation::Codon_translation()
{
}



void Codon_translation::define_translation_tables()
{

    string aa[167] =   {"M",
                        "W",
                        "F", "F", "F",
                        "Y", "Y", "Y",
                        "C", "C", "C",
                        "H", "H", "H",
                        "Q", "Q", "Q",
                        "N", "N", "N",
                        "K", "K", "K",
                        "D", "D", "D",
                        "E", "E", "E",
                        "I", "I", "I", "I", "I", "I", "I",
                        "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P",
                        "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
                        "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V",
                        "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
                        "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G",
                        "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S",
                        "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L",
                        "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R",
                        "X", "-"
                       };

    string cod[167] =  {"ATG",
                        "TGG",
                        "TTT", "TTC", "TTY",
                        "TAT", "TAC", "TAY",
                        "TGT", "TGC", "TGY",
                        "CAT", "CAC", "CAY",
                        "CAA", "CAG", "CAR",
                        "AAT", "AAC", "AAY",
                        "AAA", "AAG", "AAR",
                        "GAT", "GAC", "GAY",
                        "GAA", "GAG", "GAR",
                        "ATT", "ATC", "ATH", "ATA", "ATY", "ATW", "ATM",
                        "CCT", "CCC", "CCA", "CCG", "CCN", "CCY", "CCR", "CCM", "CCK", "CCS", "CCW", "CCB", "CCD", "CCH", "CCV",
                        "ACT", "ACC", "ACA", "ACG", "ACN", "ACY", "ACR", "ACM", "ACK", "ACS", "ACW", "ACB", "ACD", "ACH", "ACV",
                        "GTT", "GTC", "GTA", "GTG", "GTN", "GTY", "CTR", "GTM", "GTK", "GTS", "GTW", "GTB", "GTD", "GTH", "GTV",
                        "GCT", "GCC", "GCA", "GCG", "GCN", "GCY", "GCR", "GCM", "GCK", "GCS", "GCW", "GCB", "GCD", "GCH", "GCV",
                        "GGT", "GGC", "GGA", "GGG", "GGN", "GGY", "GGR", "GGM", "GGK", "GGS", "GGW", "GGB", "GGD", "GGH", "GGV",
                        "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "TCN", "TCY", "TCR", "TCM", "TCK", "TCS", "TCW", "TCB", "TCD", "TCH", "TCV", "AGY",
                        "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "CTN", "CTY", "CTR", "CTM", "CTK", "CTS", "CTW", "CTB", "CTD", "CTH", "CTV", "TTR",
                        "CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "CGN", "CGY", "CGR", "CGM", "CGK", "CGS", "CGW", "CGB", "CGD", "CGH", "CGV", "AGR",
                        "NNN", "---"
                       };


    for (int i=0; i<167; i++)
    {
        codon_to_aa.insert(make_pair(cod[i],aa[i]));
    }
}

string Codon_translation::gapped_DNA_to_protein(string *sequence) const
{
    string prot;

    for (unsigned int j=0; j<sequence->length(); j+=3)
    {
        string codon = sequence->substr(j,3);
        if (codon_to_aa.find(codon) == codon_to_aa.end())
        {
            sequence->replace(j,3,"NNN");
            prot += "X";
        }
        else
        {
            prot += codon_to_aa.find(codon)->second;
        }
    }

    return prot;
}

//"M",
//"ATG",

//"W",
//"TGG",

//"F", "F", "F",
//"TTT", "TTC", "TTY",

//"Y", "Y", "Y",
//"TAT", "TAC", "TAY",

//"C", "C", "C",
//"TGT", "TGC", "TGY",

//"H", "H", "H",
//"CAT", "CAC", "CAY",

//"Q", "Q", "Q",
//"CAA", "CAG", "CAR",

//"N", "N", "N",
//"AAT", "AAC", "AAY",

//"K", "K", "K",
//"AAA", "AAG", "AAR",

//"D", "D", "D",
//"GAT", "GAC", "GAY",

//"E", "E", "E",
//"GAA", "GAG", "GAR",

//"I", "I", "I", "I", "I", "I", "I",
//"ATT", "ATC", "ATH", "ATA", "ATY", "ATW", "ATM",

//"P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P",
//"CCT", "CCC", "CCA", "CCG", "CCN", "CCY", "CCR", "CCM", "CCK", "CCS", "CCW", "CCB", "CCD", "CCH", "CCV",

//"T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
//"ACT", "ACC", "ACA", "ACG", "ACN", "ACY", "ACR", "ACM", "ACK", "ACS", "ACW", "ACB", "ACD", "ACH", "ACV",

//"V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V",
//"GTT", "GTC", "GTA", "GTG", "GTN", "GTY", "CTR", "GTM", "GTK", "GTS", "GTW", "GTB", "GTD", "GTH", "GTV",

//"A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
//"GCT", "GCC", "GCA", "GCG", "GCN", "GCY", "GCR", "GCM", "GCK", "GCS", "GCW", "GCB", "GCD", "GCH", "GCV",

//"G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G",
//"GGT", "GGC", "GGA", "GGG", "GGN", "GGY", "GGR", "GGM", "GGK", "GGS", "GGW", "GGB", "GGD", "GGH", "GGV",

//"S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S",
//"TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "TCN", "TCY", "TCR", "TCM", "TCK", "TCS", "TCW", "TCB", "TCD", "TCH", "TCV", "AGY",

//"L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L",
//"TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "CTN", "CTY", "CTR", "CTM", "CTK", "CTS", "CTW", "CTB", "CTD", "CTH", "CTV", "TTR",

//"R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R", "R",
//"CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "CGN", "CGY", "CGR", "CGM", "CGK", "CGS", "CGW", "CGB", "CGD", "CGH", "CGV", "AGR",
