#include "nudb_kmer_db.h"
#include "kmer_nudb.h"
#include "fasta_parser.h"
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/math/statistics/linear_regression.hpp>

#include <stdexcept>

namespace fs = boost::filesystem;

int main(int argc, char **argv)
{
    if (argc < 3)
    {
	std::cerr << "Usage: nudb_call_kmers kmer.nudb function.index [file ...]\n";
	exit(1);
    }

    std::string db_file = argv[1];
    std::string function_index_file = argv[2];
    typedef NuDBKmerDb<8> KDB;
    KDB db(db_file);

    if (!db.exists())
    {
	std::cerr << "database does not exist\n";
	exit(1);
	// db.create();
    }
    db.open();

    KmerNudb kmer_nudb(db, function_index_file);

    try {
	FastaParser parser;
    
	parser.set_callback([&kmer_nudb](const std::string &id, const std::string &seq) {
	
	    double slen = static_cast<double>(seq.length());
	    auto calls = std::make_shared<std::vector<KmerCall>>();
	    // auto hits = std::make_shared<std::vector<hit_in_sequence_t>>();
	    auto ostats = std::make_shared<KmerOtuStats>();
	
	    auto hit_cb = [ slen, &kmer_nudb](const hit_in_sequence_t &hit) {
		auto kd = hit.kdata;
		if (true)
		{
		    std::cout << hit.kmer << "\t" << kmer_nudb.function_at_index(kd->function_index) << "\t" << kd->median << "\t" << kd->mean << "\t" << kd->var << "\t" << sqrt(kd->var) << "\t" << kd->n_proteins << "\n";
		}
	    };
	    kmer_nudb.process_aa_seq(id, seq, calls, hit_cb, ostats);
	    for (auto c: *calls)
	    {
		// std::cout << c << "\n";
	    }
	    FunctionIndex fi;
	    std::string func;
	    float score;
	    float wt;
	    float offset;
	    kmer_nudb.find_best_call(*calls, fi, func, score, wt, offset);
	    std::cout << id << "\t" << func << "\t" << fi << "\t" << score << "\n";
	    return 0;
	});

	if (argc == 3)
	{
	    parser.parse(std::cin);
	    parser.parse_complete();
	}
	else
	{
	    for (int i = 3; i < argc; i ++)
	    {
		fs::path p(argv[i]);
		fs::ifstream ifstr(p);
		parser.parse(ifstr);
		parser.parse_complete();
	    }
	}
    }
    catch (std::runtime_error &x)
    {
	std::cerr<< "caught " << x.what() << "\n";
    }
}

