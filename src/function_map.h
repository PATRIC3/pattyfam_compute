#ifndef _function_map_h
#define _function_map_h

/*!
 * Manage the ID to function mapping.
 *
 * FunctionMap also maintains the database of function => genome mappings
 * that is used to determine whether a given function has enough evidence
 * available to be used to create signature kmers.
 *
 * Functions occurring in genomes are determined by the presence of a
 * protein with the given function occurring in the fasta file for
 * a genome.
 *
 * The genome for a fasta file is determined using either the identifier
 * from the first sequence in the file (if it is a fig| identifier), or
 * the genome name in square brackets at the end of the definition line
 * in the case of fasta files loaded in standard genbank format from
 * outside the SEED environment.
 *
 * 		
 */
class FunctionMap
{
public:
    FunctionMap(const std::string &kept_file = "") {
	if (!kept_file.empty())
	{
	    std::cerr << "writing kept data to " << kept_file << "\n";
	    kept_function_stream_.open(kept_file);
	}
    }

    /*! @brief Load ID assignments from the given file.

       We assume the file contains tab delimited id, definition pairs .

       The assignments we store have # comments stripped.
    */
    void load_id_assignments(const fs::path &file) {
	fs::ifstream ifstr(file);
	std::string line;
	int lineno = 0;
	while (std::getline(ifstr, line))
	{
	    lineno++;
	    size_t s = line.find('\t');
	    if (s == std::string::npos)
	    {
		std::cerr << "bad line " << lineno << " in file " << file << "\n";
		continue;
	    }
	    size_t s2 = line.find('\t', s+1);
	    std::string id = line.substr(0, s);
	    std::string func;
	    if (s2 == std::string::npos)
		func = line.substr(s+1);
	    else
		func = line.substr(s+1, s2 - s - 1);

	    std::string stripped, delim, comment;
	    
	    split_func_comment(func, stripped, delim, comment);

	    original_assignment_stripped_[id] = stripped;
	    original_assignment_[id] = func;
	    //
	    // If we have an assignment marking the protein as truncated,
	    // we need to remove any existing assignment. Otherwise the protein
	    // will be used with the original assignment.
	    //
	    if (delim == "#" && is_truncated_comment(comment))
	    {
		// std::cerr << "skipping truncation " << func << "\n";
//		id_function_map_.erase(id);
		continue;
	    }
	    // std::string stripped = strip_func_comment(func);
	    id_function_map_[id] = stripped;
	    // std::cerr << "Load '" << id << "' as '" << func << "' stripped='" << stripped << "'\n";
	}
    }

    /*! @brief Load assignments and sequence visibility data from a fasta file.
     
      If the def line has an assignment on it, update the assignment map.
      
      We also update the function_genome map to record which functions appear
      in which genomes.
      
      We make the simplifying assumption that there is a 1:1 mapping
      between the proteins in a genome and the contents of a given fasta file.

      @param file Fasta file of protein data
      @param keep_function_flag If true, mark each function found in this file as a kept function.
      @param deleted_fids List of protein identifiers to exclude from the mapping.
     */
    void load_fasta_file(const fs::path &file, bool keep_function_flag, const std::set<std::string> &deleted_fids) {

	const boost::regex genome_regex("\\s+(.*)\\s+\\[([^]]+)\\]$");
	const boost::regex figid_regex("fig\\|(\\d+\\.\\d+)");
	const boost::regex genome_id_regex("\\d+\\.\\d+");
	
	fs::ifstream ifstr(file);

	FastaParser parser;

	std::string genome;

	parser.set_def_callback([this, &deleted_fids, &genome, &genome_regex, &figid_regex, &genome_id_regex, &file, keep_function_flag]
				(const std::string &id, const std::string &def, const std::string &seq) {
		if (id.empty())
		    return;
		else if (deleted_fids.find(id) != deleted_fids.end())
		    return;

		boost::smatch match;

		//
		// Need to always parse for [genome]
		//

		std::string func;
		if (!def.empty())
		{
		    size_t x = def.find_first_not_of(" \t");
		    func = def.substr(x);
		}
		std::string genome_loc;
		if (boost::regex_match(def, match, genome_regex))
		{
		    std::string delim, comment;
		    split_func_comment(match[1], func, delim, comment);
		    if (delim == "#" && is_truncated_comment(comment))
		    {
			// std::cerr << "skipping truncation " << match[1] << "\n";
//			id_function_map_.erase(id);
			return;
		    }
		    //func = strip_func_comment(match[1]);
		    genome_loc = match[2];
		}
		
		// Determine genome from first sequence.
		if (genome.empty())
		{
		    if (def.empty())
		    {
			if (boost::regex_search(id, match, figid_regex))
			{
			    genome = match[1];
			}
		    }
		    else
		    {
			if (!genome_loc.empty())
			{
			    genome = genome_loc;
			}
		    }
		}
		if (genome.empty())
		{
		    // default it to the file, just to have a value
		    genome = file.filename().string();
		    
		    if (!boost::regex_match(genome, genome_id_regex))
		    {
			std::cerr << "cannot determine genome from file " << file << "\n";
		    }
		}

		/*
		 *
		 * If the current-function map has a value, use that. We assume explicit setting
		 * of functions overrides what is in the fasta.
		 *
		 * If we're assigning a function, update the id to function
		 * map.
		 *
		 * Then look up the function for this id and add to the function_genome map.
		 */
		std::string cur_func = id_function_map_[id];
		if (cur_func.empty())
		{
		    if (!func.empty())
		    {
			id_function_map_[id] = func;
		    }
		}
		else
		{
		    func = cur_func;
		}
		
		if (func.empty())
		{
		    // std::cerr << "No function found for " << id << "\n";
		}
		else
		{
		    function_genome_map_[func].insert(genome);
		    if (keep_function_flag)
		    {
			// std::cerr << "Keeping function " << func << "\n";
			good_functions_.insert(func);
		    }
		}
			
		return;
	    });
	parser.parse(ifstr);
	parser.parse_complete();
    }

    /*!
     * @brief Process the contents of function_genome_map to determine which
     * function for which we should build signatures.
     *
     * Criteria - one of the following must be met:
     *
     * *   Genome count for a function >= min_reps_required
     *
     * *   One of the roles in the function is in the good_roles set.
     *
     * *   The function is in the good_functions set.
     *
     * For each function that we keep we will assign an identifier in the
     * function_index_map which is then used as the means to describe
     * which functions are to have signatures created.
     * @param min_reps_required Function must occur in this many distinct genomes to be kept
     */
    void process_kept_functions(int min_reps_required, const std::set<std::string> &ignored_functions) {
	std::set<std::string> kept;
	for (auto entry: function_genome_map_)
	{
	    auto function = entry.first;
	    int n_genomes = (int) entry.second.size();
	    if (kept_function_stream_.is_open())
		kept_function_stream_ << function << ": " << n_genomes << " genomes\n";
	    bool ok = false;

	    if (n_genomes >= min_reps_required)
	    {
		if (kept_function_stream_.is_open())
		    kept_function_stream_ << "Keeping " << function << ": enough genomes\n";
		ok = true;
	    }
	    else if (good_functions_.find(function) != good_functions_.end())
	    {
		if (kept_function_stream_.is_open())
		    kept_function_stream_ << "Keeping " << function << ": in good functions list\n";
		ok = true;
	    }
	    else
	    {
		std::vector<std::string> roles = roles_of_function(function);
		if (kept_function_stream_.is_open())
		    kept_function_stream_ << "Role check " << function << ":\n";
		for (auto role: roles)
		{
		    if (good_roles_.find(role) != good_roles_.end())
		    {
			if (kept_function_stream_.is_open())
			    kept_function_stream_ << "  Keeping " << function << ": " << role << " in good roles list\n";
			ok = true;
			break;
		    }
		    else
		    {
			if (kept_function_stream_.is_open())
			    kept_function_stream_ << "  " << function << ": " << role << " not in list\n";
		    }
		}

		if (!ok)
		{
		    if (kept_function_stream_.is_open())
			kept_function_stream_ << "Reject " << function << "\n";
		}
	    }
	    if (ok)
		kept.insert(function);
	}

	/*
	 * Remove functions that we are going to ignore.
	 */
	for (auto fn: ignored_functions)
	{
	    std::cerr << "Ignore '" << fn << "'\n";
	    kept.erase(fn);
	}

	/*
	 * Assign sequential function IDs.
	 */
	unsigned short next = 0;
	for (auto f: kept)
	{
	    unsigned short id = next++;
	    function_index_map_[f] = id;
	    index_function_map_[id] = f;
	}
	std::cout << "kept " << next << " functions\n";
    }

    void dump() {
	std::ofstream ofstr("fm.dump");
	ofstr << "function_genome_map\n";
	for (auto x: function_genome_map_)
	{
	    ofstr << x.first << ":";
	    for (auto y: x.second)
		ofstr << " " << y;
	    ofstr << "\n";
	}
	ofstr << "id_function_map\n";
	for (auto x: id_function_map_)
	{
	    ofstr << x.first << " '" << x.second << "'\n";
	}
    }

    void lookup_original_assignment(const std::string &id, std::string &func, std::string &stripped) {
	func = original_assignment_[id];
	stripped = original_assignment_stripped_[id];
    }

    std::string lookup_function(FunctionIndex idx) {
	auto it = index_function_map_.find(idx);
	if (it == index_function_map_.end())
	    return "";
	else
	    return it->second;
    }

    std::string lookup_function(const std::string &id) {
	auto it = id_function_map_.find(id);
	if (it == id_function_map_.end())
	    return "";
	else
	    return it->second;
    }
    FunctionIndex lookup_index(const std::string &func) {
	auto it = function_index_map_.find(func);
	if (it == function_index_map_.end())
	    return USHRT_MAX;
	else
	    return it->second;
    }

    void write_function_index(const fs::path &dir) {
	fs::ofstream of(dir / "function.index");

	std::map<int, std::string> by_index;
	for (auto ent: function_index_map_)
	    by_index.insert(std::make_pair(ent.second, ent.first));
	for (auto ent: by_index)
	{
	    of << ent.first << "\t" << ent.second << "\n";
	}
    }

    void add_good_roles(const std::vector<std::string> &r) {
	std::copy(r.begin(), r.end(), std::inserter(good_roles_, good_roles_.end()));
    }
    
    void add_good_functions(const std::vector<std::string> &r) {
	std::copy(r.begin(), r.end(), std::inserter(good_functions_, good_functions_.end()));
    }
    
private:
    /*!
      @brief Mapping from function string to the list of genomes in which it occurs.

      Initialized from scanning the protein fasta files, which we assume
      are named with the genome ID they came from.
    */
    std::map<std::string, std::set<std::string> > function_genome_map_;

    /*! @brief Mapping from protein identifier to assigned function.

      Initialized from the function definition files; assignments
      are overridden if a fasta file also has a function definition in the
      fasta header.
    */
    std::map<std::string, std::string> id_function_map_;
    std::map<std::string, FunctionIndex> function_index_map_;
    std::map<FunctionIndex, std::string> index_function_map_;

    std::set<std::string> good_roles_;
    std::set<std::string> good_functions_;

    std::ofstream kept_function_stream_;

    std::map<std::string, std::string> original_assignment_stripped_;
    std::map<std::string, std::string> original_assignment_;
};


#endif
