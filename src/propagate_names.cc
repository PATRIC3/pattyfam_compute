
#include "propagate_names.h"
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"

#include <thread>
#include <vector>
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>

#include <unistd.h>

// include log4cxx header files.
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "log4cxx/propertyconfigurator.h"
#include "log4cxx/helpers/exception.h"

#include "operators.h"

using namespace log4cxx;
using namespace log4cxx::helpers;
// Define a static logger variable so that it references the
// Logger instance named "MyApp".
LoggerPtr logger(Logger::getLogger("propagate_names"));

using namespace boost::filesystem;
namespace po = boost::program_options;

void FamData::read_pegsyn_file(path &data_file)
{
    std::ifstream ifile(data_file.string());
    std::string line;

    while (std::getline(ifile, line))
    {
	//std::vector<std::string> cols;
	//boost::split(cols, line, boost::is_any_of("\t"));
	// gnl|md5|8eea9950bc1e0f99f75a4a976304d636,169	fig|438.15.peg.2419,169;fig|940265.3.peg.655,169;fig|940265.4.peg.1650,169

	if (line.substr(0, 8) != "gnl|md5|")
	    throw std::runtime_error("Invalid pegsyn line");

	size_t com = line.find(",", 8);
	size_t tab = line.find("\t", com + 1);

	if (com == std::string::npos)
	    throw std::runtime_error("Invalid pegsyn line (no comma)");

	std::string md5(line.substr(8, com - 8));

	std::string rest(line.substr(tab + 1));

 	// std::cout << "have md5 '" << md5 << "' \t" << rest << "\n";

	size_t pos = 0;
	bool first = true;
	while (pos < rest.size())
	{
	    size_t nxt = rest.find(",", pos);
	    if (nxt == std::string::npos)
		break;
	    std::string fid(rest.substr(pos, nxt - pos));
	    // std::cout << fid << "\n";

	    if (first)
	    {
		/*
		 * Attempt to insert the key fid. This may fail if another
		 * thread attempts to insert a key for the same md5.
		 * That is OK, since we will end up with a unique key
		 * for that md5.
		 * If we succeed, insert a record into the fid_is_key
		 * map.
		 */
		 auto result = md5_to_key_.insert(std::make_pair(md5, fid));
		 // Hm. Previous version of code always inserted to fid_is_key even if insert failed.
		 if (result.second)
		 {
		     // std::cout << "inserted key " << fid << " for " << md5 << "\n";
		     fid_is_key_.insert(std::make_pair(fid, md5));
		 }
		 else
		 {
		     LOG4CXX_DEBUG(logger, "Already key " << result.first->second << " for md5 " << md5);
		 }
		first = false;
	    }
	    fid_to_md5_.insert(std::make_pair(fid, md5));
	    
	    nxt = rest.find(";", nxt);
	    if (nxt == std::string::npos)
		break;
	    pos = nxt + 1;
	}
    }
    std::cout << "fid_is_key " << fid_is_key_.load_factor() << " md5_to_key " << md5_to_key_.load_factor()
	      << " fid_to_md5 " << fid_to_md5_.load_factor() << "\n";
    std::cout << "fid_is_key " << fid_is_key_.size() << " md5_to_key " << md5_to_key_.size()
	      << " fid_to_md5 " << fid_to_md5_.size() << std::endl;
//    std::cout << fid_is_key_<< "\n";

}

void FamData::read_pegsyn()
{
    directory_iterator end_iter;
    tbb::concurrent_vector<path> paths;
    
    for (auto gpath : directory_iterator(data_dir_))
    {
	if (is_directory(gpath))
	{
	    std::string genus(gpath.path().filename().string());

	    if (target_genus_ == "" || genus == target_genus_)
	    {
		std::cout << genus << "\n";

		path pegsyn_file = gpath.path() / "nr/peg.synonyms";

		if (!is_regular_file(pegsyn_file))
		    throw std::runtime_error("Pegsynfile " + pegsyn_file.string() + " does not exist");

		paths.push_back(pegsyn_file);
	    }
	}
    }

    tbb::parallel_for( paths.range(),
		       [this](auto r) {
		      
			   for(auto ent =r.begin(); ent!=r.end(); ++ent)
			      {
				  path pegsyn_file = *ent;
				  std::cout << "Process path " << pegsyn_file.string() << "\n";;
				  try {
				      this->read_pegsyn_file(pegsyn_file);
				  } catch (std::runtime_error e)
				  {
				      std::cerr << "couldn't process " << pegsyn_file.string() << ": " << e.what() << "\n";
				  }
				  std::cout << "Done with " << pegsyn_file.string() << "\n";;
			      }
		       }
	);
} 
   

void FamData::read_fams_file()
{
    std::ifstream ifile(fams_file_.string());
    std::string line;

    std::string last_fam;

    int line_number = 0;
    while (std::getline(ifile, line))
    {
	if (line_number++ % 1000000 == 0)
	{
	    std::cout << fams_file_.string() << " line " << line_number
		      << " md5_to_fam_: " << md5_to_fam_.load_factor()
		      << " fam_to_md5s_: " << fam_to_md5s_.load_factor()
		      << " fam_to_function_ " << fam_to_function_.load_factor() << std::endl;
	}

	// GF00000000	50	38	fig|1787.5.peg.2249	242	(2-pyrone-4,6-)dicarboxylic acid hydrolase	53985	Mycobacterium 53985
	// 0	gfam id
	// 1	fams merged
	// 2	genera merged
	// 3	feature id
	// 4	feature length
	// 5	function
	// 6	local fam
	// 7	genus
	// 8	local fam (redundant)

	std::vector<std::string> cols;
	boost::split(cols, line, boost::is_any_of("\t"));

	std::string &peg(cols[3]);

	bool debug = cols[0] == "GF00095422";

	/*
	 * Look up md5 for this member.  We don't process more than
	 * one entry in a family with the given md5.
	 */
	    
	#if 0
	// after changing fam_to_md5s_ to a set, don't need this
	auto fkiter = fid_is_key_.find(peg);
	if (fkiter == fid_is_key_.end())
	{
	    if (debug)
		std::cerr << "Skipping " << peg << " due to non-key" << std::endl;
	    continue;
	}
	#endif

	auto fkiter = fid_to_md5_.find(peg);
	if (fkiter == fid_to_md5_.end())
	{
	    std::cerr << "Could not find md5 for " << peg << std::endl;
	    continue;
	}
	std::string &md5 = fkiter->second;
	if (debug)
	    std::cerr << "Key: " << peg << " " << md5 << std::endl;

	std::string &gf(cols[0]);
	std::string &function(cols[5]);
	std::string &lnum(cols[6]);
	std::string &genus(cols[7]);

	std::string fam;
	if (family_type_ == FAM_GLOBAL)
	{
	    fam = gf;
	}
	else
	{
	    fam = genus + "." + lnum;
	}
	if (fam != last_fam)
	{
	    fam_to_function_.insert(std::make_pair(fam, function));
	    last_fam = fam;
	}
	auto check = md5_to_fam_.insert(std::make_pair(md5, fam));
	if (!check.second)
	{
	    LOG4CXX_DEBUG(logger, "Peg " << peg << " md5 " << md5 << " in fam " << fam << " already in fam " << check.first->second);
	}
	    
	md5_to_fam_.insert(std::make_pair(md5, fam));
	//fam_to_md5s_.insert(std::make_pair(fam, md5));
	fam_to_md5s_[fam].insert(md5);
	    
//	std::cout << line << "\n";
    }
    
}

RenumberState::RenumberState(FamData &old_data, FamData &new_data)
    : old_data_(old_data), new_data_(new_data), new_idx_(1)
{
}


void RenumberState::phase_1_body(const std::string &fam, std::set<std::string> &fids)
{
    LOG4CXX_DEBUG(logger, "Phase 1 body for " + fam);

    std::set<std::string> nfam_checked;
    std::map<std::string, int> nfam_count;
    int bad = 0;
    
    for (auto peg: fids)
    {
	if (!new_data_.exists(peg))
	    continue;

	auto nfam = new_data_.peg_to_fam(peg);
	auto nfun = new_data_.fam_to_fun(nfam);

	/*
	 * Verify we haven't looked at the given new family yet.
	 *
	 * The size of nfam_checked is the number of new fams that
	 * this old fam maps to. Save this, because we need it in
	 * phase 2.
	 *
	 * We actually need to save the set, because the operation
	 * we use later collects a set of sets across a number of npegs.
	 */

	auto check = nfam_checked.emplace(nfam);
	if (!check.second)
	    continue;

	LOG4CXX_DEBUG(logger, peg << ": " << nfam << " " << nfun);
	
	for (auto npeg: new_data_.fam_to_md5s_[nfam])
	{
	    if (old_data_.exists(npeg))
	    {
		auto npeg_ofam = old_data_.peg_to_fam(npeg);
		if (npeg_ofam == fam)
		{
		    nfam_count[nfam]++;
		}
		else
		{
		    bad++;
		    if (bad > 10)
		    {
			LOG4CXX_DEBUG(logger, "exiting for bad=" << bad);
			break;
		    }
		}	
	    }
	}
    }
    old_fam_to_new_fam_set_.emplace(fam, nfam_checked);
    for (auto nf: nfam_count)
    {
	LOG4CXX_DEBUG(logger, "  " << nf.first << " " << nf.second);
    }

    if (!bad)
    {
	if (nfam_count.size() == 1)
	{
	    auto nfam = nfam_count.begin()->first;
	    log_result(nfam + " NOW " + fam + "\n");

	    new_fam_name_[nfam] = fam;
	    old_fam_used_[fam] = nfam;
	}
	else if (nfam_count.size() > 1)
	{
	    /*
	      std::cerr << "split " << fam << ":\n";
	    for (auto x: nfam_count)
	    {
		std::cerr << "  " << x.first << " " << x.second << "\n";
	    }
	    */

	    // Sort fams by count
	    std::vector<std::pair<std::string, int>> vec = sort_by_values(nfam_count);
	    std::stringstream ss;
	    ss << "SPLIT O " << fam << " => N";
	    for (auto x: vec)
		ss << " " << x.first;
	    ss << "\n";

	    log_result(ss.str());

	    /*
	     * On a split, we propagate the first name and allocate
	     * new names for the rest.
	     */

	    auto fiter = vec.begin();
	    auto nfam = fiter->first;
	    new_fam_name_[nfam] = fam;
	    old_fam_used_[fam] = nfam;
	    fiter++;

	    std::string res = nfam + " NOW " + fam + "\n";
	    log_result(res);

	    while (fiter != vec.end())
	    {
		std::string nm = allocate_new_id();
		new_fam_name_[fiter->first] = nm;
		std::string res = fiter->first + " NOW " + nm + "\n";
		log_result(res);
		fiter++;
	    }
	}
    }
}

void RenumberState::phase_1()
{
    LOG4CXX_INFO(logger, "Begin phase 1 fam count " << old_data_.fam_to_md5s_.size());

    if (1)
    {
	tbb::parallel_for(old_data_.fam_to_md5s_.range(),
			  [&](string_set_map_t::range_type &r)
			  {
			      LOG4CXX_DEBUG(logger,"Inside parallel for ");
			      
			      for (auto &elt: r)
			      {
				  phase_1_body(elt.first, elt.second);
			      }
			  }
	    );
    }
    else
    {
	for (auto &elt: old_data_.fam_to_md5s_)
	{
	    phase_1_body(elt.first, elt.second);
	}
    }
    LOG4CXX_INFO(logger, "End phase 1");
}

void RenumberState::phase_2_body(const std::string &nfam, std::set<std::string> &nfids)
{
    LOG4CXX_INFO(logger, "Phase 2 body for " << nfam);

    if (new_fam_name_.find(nfam) != new_fam_name_.end())
    {
	LOG4CXX_INFO(logger, "Phase 2 body for " << nfam << " done, no members found");
	return;
    }

    std::map<std::string, int> mapped_nfams;
    std::map<std::string, int> ocount;
    std::vector<std::string> npegs_that_exist(nfids.size());

    auto iter = std::copy_if(nfids.begin(), nfids.end(), npegs_that_exist.begin(),
			     [this](const std::string &fid) { return old_data_.exists(fid); });
    npegs_that_exist.resize(std::distance(npegs_that_exist.begin(),iter));
    
    if (npegs_that_exist.empty())
    {
	LOG4CXX_DEBUG(logger, "  all new pegs");
	std::string nm = allocate_new_id();
	new_fam_name_[nfam] = nm;
	log_result(nfam + " NOW " + nm + "\n");
    }
    else
    {
	LOG4CXX_INFO(logger, "starting npegs_that_exist loop " << npegs_that_exist.size());

	for (auto &npeg: npegs_that_exist)
	{
	    /*
	      
	      For each peg in the new family, if the corresponding peg exists
	      in the old families then $ofam is that family.
	      
	      We examine each peg in the old family ($opeg below)
	      and map it back to the corresponding new family.
	      We maintain a count of the new new families thus mapped.
	      
	      If all of those mappings fold back to the same new family $nfam
	      then we have a case where the old families joined to form the new family.
	    */

	    auto ofamit = old_data_.md5_to_fam_.find(npeg);
	    if (ofamit == old_data_.md5_to_fam_.end())
	    {
		LOG4CXX_DEBUG(logger, "No mapping for " << npeg << " found in old fams");
	    }
	    else
	    {
		// TODO tweak this logic for making md5s a map 
		auto &ofam(ofamit->second);
		
		// LOG4CXX_DEBUG(logger, "  O " << npeg << " " << old_data_.md5_to_key_[npeg] << " " << ofam << " " << old_data_.fam_to_function_[ofam]);

		if (ocount[ofam] == 0)
		{
		    /*
		     * Look up the count of new fams that this old fam
		     * maps to, as computed in phase 1.
		     */
		    #if 1
		    auto iter = old_fam_to_new_fam_set_.find(ofam);
		    
		    if (iter != old_fam_to_new_fam_set_.end())
		    {
			for (auto mapped_nfam: iter->second)
			    mapped_nfams[mapped_nfam]++;
		    }
		    
		    #else
                    auto list_iter = old_data_.fam_to_md5s_.find(ofam);
		    if (list_iter != old_data_.fam_to_md5s_.end())
		    {
			auto &list = list_iter->second;
			LOG4CXX_INFO(logger, "process list for " << ofam << " " << list.size());
			for (auto opeg: list)
			{
			    if (!new_data_.exists(opeg))
				continue;
			    auto mapped_nfam = new_data_.peg_to_fam(opeg);
			    // LOG4CXX_DEBUG(logger, "  O " << opeg << " " << old_data_.md5_to_key_[opeg] << " N " << mapped_nfam << " " << new_data_.fam_to_function_[mapped_nfam]);
			    mapped_nfams[mapped_nfam]++;
			}
		    }
		    #endif
		}
		ocount[ofam]++;
	    }
	}
	LOG4CXX_INFO(logger, "found mapped_nfams " << mapped_nfams.size());

	if (mapped_nfams.size() == 1)
	{
	    std::vector<std::pair<std::string, int>> ocount_sorted = sort_by_values(ocount);
	    std::string rest;
	    for (auto x: ocount_sorted)
	    {
		if (rest != "")
		{
		    rest += " ";
		}
		rest += x.first;
		LOG4CXX_DEBUG(logger, x.first << ": " << x.second << " " << old_data_.fam_to_function_[x.first]);
	    }
	    auto oname = ocount_sorted[0].first;
	    
	    new_fam_name_[nfam] = oname;
	    old_fam_used_[oname] = nfam;
	    log_result(nfam + " NOW " + oname + "\n");
	    log_result("JOIN " + rest + " => " + nfam + "\n");
	}
    }
    LOG4CXX_INFO(logger, "Phase 2 body for " << nfam << " done");
}

void RenumberState::phase_2()
{
    LOG4CXX_INFO(logger, "Begin phase 2 fam count " << new_data_.fam_to_md5s_.size());

    if (1)
    {
	tbb::parallel_for(new_data_.fam_to_md5s_.range(),
			  [&](string_set_map_t::range_type &r)
			  {
			      LOG4CXX_DEBUG(logger,"Inside parallel for ");

			      for (auto iter = r.begin(); iter != r.end(); iter++)
			      {
				  size_t i =  std::distance(iter,  r.end());
				  if (i % 1000 == 0)
				  {
				      std::cout << std::this_thread::get_id() << ": "  << i << std::endl;
				  }
				  phase_2_body(iter->first, iter->second);
			      }
			  }
	    );
    }
    else
    {
	for (auto &elt: new_data_.fam_to_md5s_)
	{
	    phase_2_body(elt.first, elt.second);
	}
    }
    LOG4CXX_INFO(logger, "End phase 2");
}

void RenumberState::phase_3_body(const std::string &fam, std::set<std::string> &fids)
{
    LOG4CXX_DEBUG(logger, "Phase 3 body for " << fam);

    if (old_fam_used_.find(fam) != old_fam_used_.end())
	return;

    std::map<std::string, int> nfams;
    int n = 0;
    for (auto fid: fids)
    {
	if (!new_data_.exists(fid))
	    continue;
	auto nfam = new_data_.peg_to_fam(fid);
	nfams[nfam]++;
	n++;
    }

    if (n == 0)
    {
	LOG4CXX_DEBUG(logger, "  No pegs, skipping");
	return;
    }

    float thresh = 0.75;
    std::vector<std::pair<std::string, int>> by_weight = sort_by_values(nfams);

    for (auto nelt: by_weight)
    {
	auto &nfam = nelt.first;
	float r = (float) nelt.second / (float) n;
	auto &renamed = new_fam_name_[nfam];
	if (renamed == "")
	{
	    LOG4CXX_DEBUG(logger, "  " << r << ": " << nfam << " " << new_data_.fam_to_function_[nfam]);
	}
	else
	{
	    LOG4CXX_DEBUG(logger, "  " << r << ": " << nfam << " " << new_data_.fam_to_function_[nfam] << " already renamed: " << renamed);
	}
    }
    auto &cand = by_weight[0].first;
    float frac = (float) by_weight[0].second / (float) n;
    auto &renamed = new_fam_name_[cand];
    if (frac > thresh && renamed == "")
    {
	new_fam_name_[cand] = fam;
	old_fam_used_[fam] = cand;
	std::stringstream ss;
	ss << cand << " NOW " << fam << " weight=" << frac << "\n";
	log_result(ss.str());
    }
}

void RenumberState::phase_3()
{
    LOG4CXX_INFO(logger, "Begin phase 3 fam count " << old_data_.fam_to_md5s_.size());

    /*
     * It is unclear whether it is safe to run this phase in parallel, since
     * one iteration may set a value that another iteration depends on.
     * So run sequentially. This phase is typically small anyway.
     */

    for (auto elt: old_data_.fam_to_md5s_)
    {
	phase_3_body(elt.first, elt.second);
    }
    LOG4CXX_INFO(logger, "End phase 3");
}

void RenumberState::write_unmapped()
{
    /*
     * Emit all of the new families and functions that have not
     * been mapped from the old families.
     *
     * Unmapped new:
     * \t fam \t function
     */

    log_result("Unmapped new:\n");
    for (auto ent: new_data_.fam_to_md5s_)
    {
	const std::string &new_fam = ent.first;

	auto fn_iter = new_fam_name_.find(new_fam);
	if (fn_iter == new_fam_name_.end() || fn_iter->second.empty())
	{	
	    log_result("\t" + new_fam + "\t" + new_data_.fam_to_function_[new_fam] + "\n");
	}
	else
	{
	    log_result("M\t" + new_fam + "\t" + new_data_.fam_to_function_[new_fam] + "\t" + fn_iter->second + "\n");
	}
	
    }

    
}

void RenumberState::set_logfile(const std::string &file)
{
    log_.open(file);
    if (log_.fail())
    {
	throw std::runtime_error("error opening " + file + ": " + strerror(errno));
    }
}

void RenumberState::logger_main()
{
    while (1)
    {
	std::string item;
	result_queue_.pop(item);
	if (item.empty())
	{
	    std::cerr << "shutting down logger\n";
	    break;
	}
	log_ << item;
    }
}

void RenumberState::start_logger()
{
    std::thread logger_thread(&RenumberState::logger_main, this);
    result_logger_thread_ = std::move(logger_thread);
}

void RenumberState::stop_logger()
{
    result_queue_.push("");
    result_logger_thread_.join();
}

void RenumberState::log_result(const std::string &res)
{
    //std::cout << res;
    //results_.push_back(res);
    result_queue_.push(res);
}
    
std::string RenumberState::allocate_new_id()
{
    std::stringstream ss;
    ss << "NEW_" << new_idx_;
    new_idx_++;
    return ss.str();
}

int main(int argc, char* argv[])
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options] fam-type old-fams old-data new-fams new-data\nAllowed options";
    po::options_description desc(x.str());

    int n_threads;
    bool debug_level;
    bool info_level;
    std::string genus;
    std::string fam_type;
    std::string old_fams, old_data;
    std::string new_fams, new_data;
    std::string log_file;
    
    desc.add_options()
	("help,h", "show this help message")
	("n-threads", po::value<int>(&n_threads)->default_value(1), "number of processing threads")
	("log-file", po::value<std::string>(&log_file), "log output to this file")
	("genus", po::value<std::string>(&genus), "genus to process")
	("debug", po::bool_switch(&debug_level), "enable debugging output")
	("info", po::bool_switch(&info_level), "enable info output")
	("fam-type", po::value<std::string>(&fam_type)->required(), "family type to process (local or global)")
	("old-fams", po::value<std::string>(&old_fams)->required(), "old family file")
	("old-data", po::value<std::string>(&old_data)->required(), "old family genus.data directory")
	("new-fams", po::value<std::string>(&new_fams)->required(), "new family file")
	("new-data", po::value<std::string>(&new_data)->required(), "new family genus.data directory")
	;
    po::positional_options_description pd;
    pd.add("fam-type", 1)
	.add("old-fams", 1)
	.add("old-data", 1)
	.add("new-fams", 1)
	.add("new-data", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(pd).run(), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
	return 1;
    }

    try {
	po::notify(vm);
    } catch (po::required_option &e)
    {
	std::cerr << "Invalid command line: " << e.what() << "\n";
	std::cerr << desc << "\n";
	exit(1);
    }

    BasicConfigurator::configure();

    logger->setLevel(Level::getWarn());
    
    if (debug_level)
	logger->setLevel(Level::getDebug());
    if (info_level)
	logger->setLevel(Level::getInfo());
    
    tbb::task_scheduler_init sched_init(n_threads);

    FamData::FamilyType ftype;
    if (fam_type == "local")
	ftype = FamData::FAM_LOCAL;
    else if (fam_type == "global")
	ftype = FamData::FAM_GLOBAL;
    else
    {
	std::cerr  << "Invalid family type (should be local or global)\n";
	exit(0);
    }
    
    path old_data_path(old_data);
    path old_fams_file(old_fams);

    if (!is_regular_file(old_fams_file))
    {
	LOG4CXX_FATAL(logger, "Cannot read old fams file " << old_fams_file);
	exit(1);
    }
    FamData old_fam_data(old_fams_file, old_data_path, genus, ftype);

    std::thread t1([&]() {
	    old_fam_data.read_pegsyn();
	    old_fam_data.read_fams_file();
	});

    path new_data_path(new_data);
    path new_fams_file(new_fams);
    FamData new_fam_data(new_fams_file, new_data_path, genus, ftype);
    std::thread t2([&]() {
	    new_fam_data.read_pegsyn();
	    new_fam_data.read_fams_file();
	});
    t1.join();
    t2.join();

    RenumberState renumber(old_fam_data, new_fam_data);

    if (log_file != "")
    {
	renumber.set_logfile(log_file);
	renumber.start_logger();
    }

    // std::cout << "Old data:\n" << old_fam_data;
    // std::cout << "New data:\n" << new_fam_data;

    
    renumber.phase_1();
    renumber.phase_2();
    renumber.phase_3();
    renumber.write_unmapped();

    if (log_file != "")
	renumber.stop_logger();
	
}
