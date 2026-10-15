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
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <iostream>
#include <sstream>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <dirent.h>
#ifdef __linux__
#include <sys/prctl.h>
#endif
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


/*****************************************************************************/
/*  --independent-queries                                                    */
/*****************************************************************************/
//
// By default the queries of --queryfile are placed one after another into
// ONE growing tree and alignment, so where a query goes and how it is aligned
// depends on the queries before it in the same file. Placing a query
//
//   - inserts a new internal node "#n#" (the target node and the query)
//     into the tree, and re-points the parent's alignment sites to it
//     (Reads_aligner::correct_sites_index()), so the next query is aligned
//     against a grown alignment;
//   - re-reconstructs the parent's ancestral sequence
//     (reconstruct_one_parsimony_ancestor()), which later queries are
//     scored and aligned against;
//   - with target preselection, adds the query and the new node to the
//     candidate targets of every later query (target_sequences);
//   - consumes an internal node number, so later nodes are named differently.
//
// With --independent-queries the reference alignment, the tree and the model
// are set up ONCE, and each query is then placed and aligned in a child
// process forked from that state: the child aligns that query alone, exactly
// as a run given only that query would, and writes what such a run writes,
// under the output name "<outfile>.<n>" (n = the query's 1-based number among
// the records READ from --queryfile: a record without a sequence is dropped by
// the reader and gets no number, which is warned about once; a repeated name
// is renamed "<name>.1", ... as in a default run). The fork's copy of the
// process is what gives every query the pristine reference: nothing a child
// changes reaches the parent or the next child.
//
// THE MANIFEST, "<outfile>.queries.tsv": a header line, one row per query in
// file order, and -- only once every query has a row -- a trailer line:
//
//   query  name  status  detail  alignment
//   ...
//   # complete: N queries, P placed, U unplaced, F failed
//
// The trailer has no tab, so a reader of five tab-separated columns skips it;
// a manifest without it was cut short (the parent was killed). In the name
// and detail columns a tab, CR or LF is written as a space (the reader splits
// a FASTA header on spaces only, so a name can hold a tab). status is
//
//   placed    the child placed the query and the alignment file named in the
//             last column exists and is not empty;
//   unplaced  the query matched no node or failed the overlap/identity
//             thresholds; no file, as a single-query run writes none;
//   failed    the child exited non-zero, was killed, or said "placed" without
//             leaving a non-empty alignment file; or the query was not run
//             (see below). detail says which.
//
// The child's outputs are the regular files "<outfile>.<n>.*" and
// "<outfile>.<n>_*". Before a query runs, such files left by an earlier run
// are removed (so a file cannot pass for this run's); after a failed query,
// what the child left is removed, and anything that cannot be is named in
// detail. A failed query loses that query only: the rest are still aligned,
// and the run exits 1. Read the manifest rather than listing the directory:
// files of an earlier, longer run under the same -o (n beyond this run's
// last query) are not touched.
//
// If pipe() or fork() fails (no memory, no processes left), nothing later
// would succeed either: that query and every later one are "failed" without
// another attempt, and the run exits 1.
//
// The parent must be single-threaded when it forks, and is: with a reference
// alignment nothing before query placement starts a thread.
//
// rand() names the scratch files of the programs a query runs (Exonerate's
// q<r>.fas, t<r>.fas). A child keeps the state main() seeded for the process:
// the children of one run, which run one after another, draw the same names
// and each removes its files before the next starts, while runs started
// apart draw different ones -- exactly as separate processes would.

namespace {

// Everything buffered is written now, or a child would write it again.
void flush_all_output()
{
    Log_output::flush();
    cout.flush();
    cerr.flush();
    fflush(NULL);
}

bool set_close_on_exec(int fd)
{
    int flags = fcntl(fd, F_GETFD);
    return flags != -1 && fcntl(fd, F_SETFD, flags | FD_CLOEXEC) != -1;
}

// A manifest field: a tab, CR or LF would add a column or a row.
string manifest_field(const string &s)
{
    string out = s;
    for(size_t k=0; k<out.size(); k++)
    {
        if(out[k] == '\t' || out[k] == '\r' || out[k] == '\n')
            out[k] = ' ';
    }
    return out;
}

// The directory and file name of an output path, and the way back to a path.
void split_output_path(const string &path, string *dir, string *base, string *prefix)
{
    size_t slash = path.rfind('/');
    if(slash == string::npos)
    {
        *dir = ".";
        *base = path;
        *prefix = "";
    }
    else
    {
        *dir = slash == 0 ? "/" : path.substr(0, slash);
        *base = path.substr(slash+1);
        *prefix = path.substr(0, slash+1);
    }
}

// Removes the given paths if they are regular files; returns those it could
// not remove, and anything else (a directory, a link, a device: not written
// by pagan2, so not removed), with the reason.
vector<string> remove_regular_files(const vector<string> &paths)
{
    vector<string> left;
    for(size_t k=0; k<paths.size(); k++)
    {
        struct stat st;
        if(lstat(paths[k].c_str(), &st) != 0)
            continue;
        if(!S_ISREG(st.st_mode))
            left.push_back(paths[k]+" (not a regular file)");
        else if(unlink(paths[k].c_str()) != 0)
            left.push_back(paths[k]+" ("+strerror(errno)+")");
    }
    return left;
}

// The paths "<outfile>.<n>.*" and "<outfile>.<n>_*" in the output directory,
// by n (no leading zero, as pagan2 writes it); n_only > 0 lists that n alone.
map<int, vector<string> > query_output_paths(const string &outfile, int n_only)
{
    map<int, vector<string> > found;
    string dir, base, prefix;
    split_output_path(outfile, &dir, &base, &prefix);
    base += ".";

    DIR *d = opendir(dir.c_str());
    if(d == 0)
        return found;
    struct dirent *e;
    while((e = readdir(d)) != 0)
    {
        string name = e->d_name;
        if(name.size() <= base.size() || name.compare(0, base.size(), base) != 0)
            continue;
        size_t k = base.size();
        if(name[k] < '1' || name[k] > '9')
            continue;
        int n = 0;
        while(k < name.size() && name[k] >= '0' && name[k] <= '9' && n < 100000000)
            n = n*10 + (name[k++]-'0');
        if(k >= name.size() || (name[k] != '.' && name[k] != '_'))
            continue;
        if(n_only > 0 && n != n_only)
            continue;
        found[n].push_back(prefix+name);
    }
    closedir(d);
    return found;
}

// A regular file, not empty (a link to one counts).
bool is_nonempty_file(const string &path)
{
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0;
}

// The child's whole life: align queries->at(i) alone, write its output under
// query_outfile, tell the parent "placed" or "unplaced" through status_fd --
// and, on a second line, that it could not arrange to die with the parent.
//
// It ends with _exit(), which skips the atexit handlers and static
// destructors the parent's copy of the process would otherwise run a second
// time. The paths that leave by exit() -- a fatal error inside placement or
// output -- do run them: none removes a file today, and none may.
void align_one_query_alone(Node *root, Model_factory *mf, int count,
                           Input_output_parser *iop, Fasta_reader *fr,
                           vector<Fasta_entry> *sequences,
                           vector<Fasta_entry> *queries, int i,
                           const string &query_outfile, int status_fd,
                           const string &pdeathsig_error)
{
    Settings_handle::st.set_string("outfile", query_outfile);

    vector<Fasta_entry> query(1, queries->at(i));

    Reads_aligner ra;
    ra.align(root,mf,count,&query);
    Node *query_root = ra.get_global_root();

    bool placed = query_root->get_number_of_query_leaves() > 0;
    iop->output_aligned_sequences(fr,sequences,query_root);

    flush_all_output();

    string status = placed ? "placed" : "unplaced";
    if(pdeathsig_error != "")
        status += "\nno-pdeathsig: "+pdeathsig_error;
    const char *p = status.c_str();
    size_t left = status.size();
    while(left > 0)
    {
        ssize_t n = write(status_fd, p, left);
        if(n < 0 && errno == EINTR)
            continue;
        if(n <= 0)
            _exit(1);
        p += n;
        left -= n;
    }
    _exit(0);
}

// Runs one query in a child; returns "placed", "unplaced" or "failed", with
// *detail saying why for "failed". *cannot_fork is set when pipe() or fork()
// failed, which no later query would escape either. *pdeathsig_error is the
// child's report that it could not ask to die with the parent.
string run_query_in_child(Node *root, Model_factory *mf, int count,
                          Input_output_parser *iop, Fasta_reader *fr,
                          vector<Fasta_entry> *sequences,
                          vector<Fasta_entry> *queries, int i,
                          const string &query_outfile, string *detail,
                          bool *cannot_fork, string *pdeathsig_error)
{
    int fds[2];
    if(pipe(fds) != 0)
    {
        *detail = string("pipe: ")+strerror(errno);
        *cannot_fork = true;
        return "failed";
    }
    // Not inherited by the helpers (Exonerate, MAFFT, ...) the child starts,
    // so the parent sees end-of-file when the child itself is gone.
    if(!set_close_on_exec(fds[0]) || !set_close_on_exec(fds[1]))
    {
        *detail = string("fcntl: ")+strerror(errno);
        close(fds[0]);
        close(fds[1]);
        *cannot_fork = true;
        return "failed";
    }

    flush_all_output();

    pid_t parent_pid = getpid();
    pid_t pid = fork();
    if(pid < 0)
    {
        *detail = string("fork: ")+strerror(errno);
        close(fds[0]);
        close(fds[1]);
        *cannot_fork = true;
        return "failed";
    }
    if(pid == 0)
    {
        string no_pdeathsig = "";
#ifdef __linux__
        // DIE WITH THE PARENT. A caller giving up on a run usually kills only
        // the pid it started (a timeout, a batch system, an operator) and
        // then removes the run's output and --temp-folder; a child left
        // running would compute its query for nobody, into a directory being
        // deleted. The kernel sends SIGKILL when the parent dies -- unless it
        // already has, between fork() and prctl(), which getppid() catches.
        // (A program the child starts, Exonerate or MAFFT, is not covered:
        // the death signal is not inherited across fork or posix_spawn. It
        // ends at its next write into the pipe its dead reader held, or with
        // the process group of the caller, if the caller kills the group --
        // which is why the child stays in the parent's process group.)
        if(prctl(PR_SET_PDEATHSIG, SIGKILL) != 0)
            no_pdeathsig = strerror(errno);
        if(getppid() != parent_pid)
            _exit(1);   // _exit(): see align_one_query_alone()
#endif
        // Other systems have no parent-death signal here (FreeBSD has
        // procctl(PROC_PDEATHSIG_CTL); macOS nothing comparable): a child
        // whose parent is killed runs until its query is aligned and
        // written, then dies of SIGPIPE reporting the outcome to the parent.
        close(fds[0]);
        align_one_query_alone(root,mf,count,iop,fr,sequences,queries,i,query_outfile,fds[1],no_pdeathsig);
        _exit(1);  // not reached
    }

    close(fds[1]);

    string message;
    char buf[256];
    while(true)
    {
        ssize_t n = read(fds[0], buf, sizeof(buf));
        if(n < 0 && errno == EINTR)
            continue;
        if(n <= 0)
            break;
        message.append(buf, n);
    }
    close(fds[0]);

    size_t eol = message.find('\n');
    if(eol != string::npos)
    {
        string note = message.substr(eol+1);
        message = message.substr(0, eol);
        if(note.compare(0, 14, "no-pdeathsig: ") == 0)
            *pdeathsig_error = note.substr(14);
    }

    int wstatus = 0;
    while(waitpid(pid, &wstatus, 0) < 0)
    {
        if(errno != EINTR)
        {
            *detail = string("waitpid: ")+strerror(errno);
            return "failed";
        }
    }

    stringstream ds;
    if(WIFSIGNALED(wstatus))
    {
        ds << "killed by signal " << WTERMSIG(wstatus) << " (" << strsignal(WTERMSIG(wstatus)) << ")";
        *detail = ds.str();
        return "failed";
    }
    if(!WIFEXITED(wstatus) || WEXITSTATUS(wstatus) != 0)
    {
        ds << "exit status " << WEXITSTATUS(wstatus);
        *detail = ds.str();
        return "failed";
    }
    if(message != "placed" && message != "unplaced")
    {
        *detail = "exited 0 without reporting an outcome";
        return "failed";
    }
    return message;
}

// Records of a FASTA file, counted by their '>' lines; -1 for anything that
// is not FASTA. The reader drops a record without a sequence, and the number
// n of every later query moves down by one.
int count_fasta_records(const string &file)
{
    ifstream in(file.c_str());
    if(!in)
        return -1;
    char c = in.get();
    while(in && (c==' ' || c=='\n'))
        c = in.get();
    if(!in || c != '>')
        return -1;
    int n = 1;
    string line;
    getline(in, line);
    while(getline(in, line))
    {
        if(!line.empty() && line[0] == '>')
            n++;
    }
    return n;
}

void write_manifest_line(FILE *manifest, const string &manifest_file, const string &line)
{
    if(fputs(line.c_str(), manifest) == EOF || fflush(manifest) != 0)
    {
        Log_output::write_out("Error: writing '"+manifest_file+"' failed: "+strerror(errno)+".\nExiting.\n\n",0);
        exit(1);
    }
}

// The parent: one child per query, in file order; the outcomes go to
// "<outfile>.queries.tsv". Returns main()'s exit status.
int align_queries_independently(Node *root, Model_factory *mf, int count,
                                Input_output_parser *iop, Fasta_reader *fr,
                                vector<Fasta_entry> *sequences)
{
    string outfile = "outfile";
    if(Settings_handle::st.is("outfile"))
        outfile = Settings_handle::st.get("outfile").as<string>();

    string format = "fasta";
    if(Settings_handle::st.is("outformat"))
        format = Settings_handle::st.get("outformat").as<string>();

    vector<Fasta_entry> queries;
    Reads_aligner::read_queries(&queries);

    int n_records = count_fasta_records(Settings_handle::st.get("queryfile").as<string>());
    if(n_records > (int)queries.size())
    {
        stringstream ws;
        ws << "\nWarning: records without a sequence in the query file: "
           << n_records-(int)queries.size() << " of " << n_records << ", dropped; the number n in '"
           << outfile << ".<n>' and in the manifest counts the " << queries.size() << " records read.\n";
        Log_output::write_out(ws.str(),0);
    }

    // FD_CLOEXEC: the programs a child starts must not hold the manifest.
    string manifest_file = outfile+".queries.tsv";
    int mfd = open(manifest_file.c_str(), O_WRONLY|O_CREAT|O_TRUNC, 0666);
    FILE *manifest = 0;
    if(mfd >= 0 && set_close_on_exec(mfd))
        manifest = fdopen(mfd, "w");
    if(manifest == 0)
    {
        int open_errno = errno;
        if(mfd >= 0)
            close(mfd);
        errno = open_errno;
        Log_output::write_out("Error: cannot write '"+manifest_file+"': "+strerror(errno)+".\nExiting.\n\n",0);
        exit(1);
    }
    write_manifest_line(manifest, manifest_file, "query\tname\tstatus\tdetail\talignment\n");

    // Files an earlier run left under the names this run writes: one listing
    // of the directory, removed query by query just before each one runs.
    map<int, vector<string> > earlier = query_output_paths(outfile, 0);

    int n_placed = 0, n_unplaced = 0, n_failed = 0;
    string cannot_fork_reason = "";
    bool pdeathsig_warned = false;

    for(int i=0; i<(int)queries.size(); i++)
    {
        int query_1based = i+1;

        stringstream qs;
        qs << outfile << "." << query_1based;
        string query_outfile = qs.str();
        string alignment_file = query_outfile+fr->get_format_suffix(format);

        string detail;
        string status;

        if(cannot_fork_reason != "")
        {
            status = "failed";
            detail = "not run: "+cannot_fork_reason;
        }
        else
        {
            // A file of an earlier run must not pass for this run's.
            if(earlier.count(query_1based))
                remove_regular_files(earlier[query_1based]);

            bool cannot_fork = false;
            string pdeathsig_error = "";
            status = run_query_in_child(root,mf,count,iop,fr,sequences,&queries,i,query_outfile,
                                        &detail,&cannot_fork,&pdeathsig_error);
            if(cannot_fork)
            {
                stringstream cs;
                cs << detail << " at query " << query_1based;
                cannot_fork_reason = cs.str();
                detail += "; no later query is run";
            }
            if(pdeathsig_error != "" && !pdeathsig_warned)
            {
                Log_output::write_out("\nWarning: prctl(PR_SET_PDEATHSIG) failed in a query process ("+
                                      pdeathsig_error+"): if this process is killed, its query processes run on.\n",0);
                pdeathsig_warned = true;
            }
            if(status == "placed" && !is_nonempty_file(alignment_file))
            {
                status = "failed";
                detail = "no alignment file written";
            }
            if(status == "failed")
            {
                vector<string> left = remove_regular_files(query_output_paths(outfile, query_1based)[query_1based]);
                for(size_t k=0; k<left.size(); k++)
                    detail += (k == 0 ? "; not removed: " : ", ")+left[k];
            }
        }

        if(status == "placed")
            n_placed++;
        else if(status == "unplaced")
            n_unplaced++;
        else
            n_failed++;

        stringstream row;
        row << query_1based << "\t" << manifest_field(queries.at(i).name) << "\t" << status << "\t"
            << manifest_field(detail) << "\t" << (status == "placed" ? alignment_file : "") << "\n";
        write_manifest_line(manifest, manifest_file, row.str());

        stringstream msg;
        msg << "\nQuery " << query_1based << "/" << queries.size() << " '" << queries.at(i).name << "': " << status;
        if(detail != "")
            msg << ", " << detail;
        if(status == "placed")
            msg << ", " << alignment_file;
        msg << "\n";
        Log_output::write_out(msg.str(),0);
    }

    // Every query has its row: say so, on a line a five-column reader skips.
    stringstream trailer;
    trailer << "# complete: " << queries.size() << " queries, " << n_placed << " placed, "
            << n_unplaced << " unplaced, " << n_failed << " failed\n";
    write_manifest_line(manifest, manifest_file, trailer.str());
    if(fclose(manifest) != 0)
    {
        Log_output::write_out("Error: writing '"+manifest_file+"' failed: "+strerror(errno)+".\nExiting.\n\n",0);
        exit(1);
    }

    struct rusage children;
    double children_cpu = 0;
    if(getrusage(RUSAGE_CHILDREN, &children) == 0)
        children_cpu = children.ru_utime.tv_sec + children.ru_utime.tv_usec/1e6
                     + children.ru_stime.tv_sec + children.ru_stime.tv_usec/1e6;

    stringstream ss;
    ss << "\nIndependent queries: " << queries.size() << " aligned one at a time; "
       << n_placed << " placed, " << n_unplaced << " unplaced, " << n_failed << " failed; "
       << "the query processes used " << children_cpu << " cpu sec.\n"
       << "Outcome per query: " << manifest_file << "\n";
    Log_output::write_out(ss.str(),0);

    return n_failed > 0 ? 1 : 0;
}

}

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
    /*  Or align each query alone against it: writes the output per query  */
    /***********************************************************************/

    int exit_code = 0;
    bool independent_queries = Settings_handle::st.is("independent-queries");

    if( independent_queries )
    {
        exit_code = align_queries_independently(root,&mf,count,&iop,&fr,&sequences);
    }


    /***********************************************************************/
    /*  If query sequences, add them to the alignment.                     */
    /***********************************************************************/

    if( Settings_handle::st.is("queryfile") && !independent_queries )
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

    if( !independent_queries )
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

    return exit_code;
}
