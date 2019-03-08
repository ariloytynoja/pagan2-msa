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

#include <string>
#include <vector>
#include <ctime>
#include "utils/settings.h"
#include "utils/settings_handle.h"
#include "utils/log_output.h"
#include "utils/input_output_parser.h"
#include "utils/fasta_reader.h"
#include "main/node.h"
#include "utils/model_factory.h"
#include "main/reads_aligner.h"

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

using namespace std;

using namespace ppa;

int main(int argc, char *argv[])
{


    /***********************************************************************/
    /*  Start the clock: requires somes OSX specific hacks                 */
    /***********************************************************************/

    clock_t t_start=clock();
    struct timespec tcpu_start, tcpu_finish;

    #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    tcpu_start.tv_sec = mts.tv_sec;
    tcpu_start.tv_nsec = mts.tv_nsec;
    #else
    clock_gettime(CLOCK_MONOTONIC, &tcpu_start);
    #endif



    /*****************************************************************************/
    /*  Read the command-line parameters; start the log system; reset the clock  */
    /*****************************************************************************/

    try
    {
        int rv = Settings_handle::st.read_command_line_arguments(argc, argv);
    }
    catch ( const boost::program_options::error& e ) {
            Settings_handle::st.info_noexit();
            Log_output::write_out("Error in command line arguments: "+string(e.what())+".\n\n",0);
            exit(1);
    }

    Log_output::open_stream();


    if(!Settings_handle::st.is("silent"))
    {
        Settings_handle::st.print_msg();
        time_t s_time;
        time( &s_time );
        Log_output::write_out("\nThe analysis started: " +string( asctime( localtime( &s_time ) ) )+"\n",0);
    }

    srand(time(0));
    clock_t analysis_start_time=clock();


    /***********************************************************************/
    /*  Threaded alignment: using maximum number by defualt                */
    /***********************************************************************/

    int n_threads = 1;
    int max_threads = boost::thread::hardware_concurrency();
    if(Settings_handle::st.is("threads"))
    {
        int nt = Settings_handle::st.get("threads").as<int > ();
        if(nt>0 && nt<=max_threads) { n_threads = nt; }
    }

    /***********************************************************************/
    stringstream ss;
    ss << "Running with "<<n_threads<<" threads.\n";
    Log_output::write_out(ss.str(),1);
    /***********************************************************************/





    /***********************************************************************/
    /*  Read the sequence file                                             */
    /***********************************************************************/

    Fasta_reader fr;
    vector<Fasta_entry> sequences;
    bool reference_alignment = false;

    Input_output_parser iop;
    iop.parse_input_sequences(&fr,&sequences,&reference_alignment);



    /***********************************************************************/
    /*  Read the guidetree file                                            */
    /***********************************************************************/

    Node *root = iop.parse_input_tree(&fr,&sequences,reference_alignment,n_threads);



    /***********************************************************************/
    /*  Check that input is fine and place the sequences to nodes          */
    /***********************************************************************/

    int data_type = -1;
    iop.match_sequences_and_tree(&fr,&sequences,root,reference_alignment,&data_type);


    /***********************************************************************/
    ss.str(string());
    ss << "Time main::input: "<<double(clock()-t_start)/CLOCKS_PER_SEC<<"\n";
    Log_output::write_out(ss.str(),"time");
    /***********************************************************************/



    /***********************************************************************/
    /*  Define the alignment model                                         */
    /***********************************************************************/

    Model_factory mf(data_type);
    iop.define_alignment_model(&fr,&mf,data_type);


    /***********************************************************************/
    ss.str(string());
    ss << "Time main::model: "<< double(clock()-t_start)/CLOCKS_PER_SEC <<"\n";
    Log_output::write_out(ss.str(),"time");
    /***********************************************************************/




    /***********************************************************************/
    /*  Read or align the sequence data.                                   */
    /***********************************************************************/

    int count = 1;
    root->name_internal_nodes(&count);

    if(reference_alignment)
    {
        root->read_reference_alignment(&mf);
    }
    else
    {
        // some features only work with 1 thread
        if(n_threads==1)
            root->start_alignment(&mf);
        else if(Settings_handle::st.is("boost"))
            root->start_threaded_alignment(&mf,n_threads);
        else
            root->start_openmp_alignment(&mf,n_threads);
    }

    /***********************************************************************/
    ss.str(string());
    ss << "Time main::align: "<< double(clock()-t_start)/CLOCKS_PER_SEC <<"\n";
    Log_output::write_out(ss.str(),"time");
    /***********************************************************************/



    /***********************************************************************/
    /*  If query sequences, add them to the alignment.                     */
    /***********************************************************************/

    if( Settings_handle::st.is("queryfile") )
    {

        Reads_aligner ra;
        ra.align(root,&mf,count);

        root = ra.get_global_root();

        /***********************************************************************/
        ss.str(string());
        ss << "Time main::reads_align: "<< double(clock()-t_start)/CLOCKS_PER_SEC <<"\n";
        Log_output::write_out(ss.str(),"time");
        /***********************************************************************/
    }



    /***********************************************************************/
    /*  Collect the results and output them                                */
    /***********************************************************************/

    iop.output_aligned_sequences(&fr,&sequences,root);


    /***********************************************************************/
    ss.str(string());
    ss << "Time main::main_exit: "<< double(clock()-t_start)/CLOCKS_PER_SEC <<"\n";
    Log_output::write_out(ss.str(),"time");
    /***********************************************************************/


    /***********************************************************************/
    /*  Stop the clock: requires somes OSX specific hacks                  */
    /***********************************************************************/

    #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    tcpu_finish.tv_sec = mts.tv_sec;
    tcpu_finish.tv_nsec = mts.tv_nsec;
    #else
    clock_gettime(CLOCK_MONOTONIC, &tcpu_finish);
    #endif

    /***********************************************************************/

    double elapsed;
    elapsed = (tcpu_finish.tv_sec - tcpu_start.tv_sec);
    elapsed += (tcpu_finish.tv_nsec - tcpu_start.tv_nsec) / 1000000000.0;

    time_t s_time;
    time( &s_time );

    /***********************************************************************/
    ss.str(string());
    ss << "\nThe analysis finished: " << asctime( localtime( &s_time ) );
    ss<<"Total time used by PAGAN: "<<elapsed<<" wall sec, " <<double_t(clock()-analysis_start_time)/CLOCKS_PER_SEC<<" cpu sec.\n\n";
    Log_output::write_out(ss.str(),0);
    /***********************************************************************/

    delete root;

    return 0;
}
