
#ifndef DEBUG_SCORING
#define DEBUG_SCORING 0
#endif


namespace fs = boost::filesystem;

// Utility for hit processing

template <class Caller>
class HitSet
{
public:
    using KData = typename Caller::KData;
    struct hit
    {
	KData kdata;
	unsigned long pos;
    };
    using HitVector = std::vector<hit>;

    HitSet(int seq_len, int min_hits, float min_weighted_hits)
	: seq_len_(seq_len)
	,min_hits_(min_hits)
	, min_weighted_hits_(min_weighted_hits)
	{
	}

    void reset();
    template< class... Args >
    void emplace_back( Args&&... args ) {
	hits_.emplace_back(std::forward<Args>(args)...);
    }

    auto rbegin() { return hits_.rbegin(); }

    bool empty() { return hits_.empty(); }
    int count() { return static_cast<int>(hits_.size()); }
    auto clear() { return hits_.clear(); }
    hit& last_hit() { return hits_.back(); }
    void process(FunctionIndex &current_fI,
		 std::shared_ptr<std::vector<KmerCall>> calls,
		 std::shared_ptr<KmerOtuStats> otu_stats) {

	int fI_count = 0;
	float weighted_hits = 0;
	typename HitVector::iterator last_hit;

	std::vector<float> protein_lengths;
	CoefficentOfDetermination codet(static_cast<double>(seq_len_));
	for (auto h_iter = hits_.begin(); h_iter != hits_.end(); h_iter++)
	{
	    if (h_iter->kdata.function_index == current_fI)
	    {
		codet.insert(static_cast<double>(h_iter->kdata.mean));
		last_hit = h_iter;
		fI_count++;
		weighted_hits += h_iter->kdata.weight;
		protein_lengths.push_back(static_cast<float>(h_iter->kdata.mean));
	    }
	}
	auto mean_length = boost::math::statistics::mean(protein_lengths);
	auto median_length = boost::math::statistics::median(protein_lengths);
	auto mad_length = boost::math::statistics::median_absolute_deviation(protein_lengths);
	// std::cout << "hit stats: " << mean_length << " " << median_length << " " << mad_length << "\n";
	// for (auto x: protein_lengths) { std::cout << x << " " ; }; std::cout << "\n";
	if ((fI_count >= min_hits_) && (weighted_hits >= min_weighted_hits_))
	{
	    if (calls)
	    {
		//R2 double r2 = codet.compute();
		// std::cerr << "computed r2=" << r2 << "\n";
		//double r2 = codet.deviation();
		calls->push_back({
			static_cast<unsigned int>(hits_[0].pos),
			static_cast<unsigned int>(last_hit->pos + (KMER_SIZE-1)),
			fI_count,
			current_fI,
			weighted_hits,
			//R2 r2,
			static_cast<unsigned int>(median_length),
			mad_length });
	    }

	    /* once we have decided to call a region, we take the kmers for fI and
	       add them to the counts maintained to assign an OTU to the sequence */
	    
	    if (otu_stats)
	    {
		for (auto ki = hits_.begin(); ki != last_hit; ki++)
		{
		    if (ki->kdata.function_index == current_fI)
		    {
			otu_stats->otu_map[ki->kdata.otu_index]++;
		    }
		}
	    }
	}
	
	auto end = hits_.rbegin();
	
	if (end[1].kdata.function_index != current_fI &&
	    end[1].kdata.function_index == end[0].kdata.function_index)
	{
	    current_fI = end[1].kdata.function_index;
	    // std::cerr << "reset cur=" << cur_fi << "\n";
	    hits_.erase(hits_.begin(), hits_.end() - 2);
	    // std::cerr << "after erase:\n";
	    //for (auto x: hits)
	    //std::cerr << x << "\n";
	}
	else {
	    hits_.clear();
	}
    }

    HitVector hits_;
    int seq_len_;
    int min_hits_;
    float min_weighted_hits_;
};

template <class Caller>
KmerGeneric<Caller>::KmerGeneric(Caller &caller, const std::string &function_index_file,
		   int min_hits, float min_weighted_hits, int max_gap) :
    caller_(caller),
    order_constraint_(false),
    min_hits_(min_hits),
    min_weighted_hits_(min_weighted_hits),
    max_gap_(max_gap)
{
    read_function_index(function_index_file);
}

template <class Caller>
KmerGeneric<Caller>::~KmerGeneric() {
}

template <class Caller>
void KmerGeneric<Caller>::read_function_index(const std::string &function_index_file)
{
    fs::ifstream ifstr(function_index_file);
    std::string line;
    int max_id = 0;
    while (std::getline(ifstr, line, '\n'))
    {
	auto tab = line.find('\t');
	int id = std::stoi(line.substr(0, tab));
	if (id > max_id)
	    max_id = id;
    }
    ifstr.close();
    ifstr.open(function_index_file);

    function_index_.resize(max_id + 1);
    
    while (std::getline(ifstr, line, '\n'))
    {
	auto tab = line.find('\t');
	int id = std::stoi(line.substr(0, tab));
	function_index_[id] = line.substr(tab + 1);
    }
}

template <class Caller>
void KmerGeneric<Caller>::process_aa_seq_hits(const std::string &id, const std::string &seq,
				   std::shared_ptr<std::vector<KmerCall>> calls,
				   std::shared_ptr<std::vector<hit_in_sequence_t<Caller>>> hits,
				   std::shared_ptr<KmerOtuStats> otu_stats)
{
    auto cb = [this, hits](hit_in_sequence_t<Caller> k) { hits->push_back(k); };
    process_aa_seq(id, seq, calls, cb, otu_stats);
}

template <class Caller>
void KmerGeneric<Caller>::gather_hits(const std::string &seqstr,
			      std::shared_ptr<std::vector<KmerCall>> calls,
			      std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
			      std::shared_ptr<KmerOtuStats> otu_stats)
{
    HitSet<Caller> hits(seqstr.length(), min_hits_, min_weighted_hits_);
    FunctionIndex current_fI = UndefinedFunction;
    double seqlen = static_cast<double>(seqstr.length());
    for_each_kmer<Caller::KmerSize>(seqstr, [this, &calls, &hit_cb, &otu_stats, &hits, &current_fI, seqlen]
				    (const std::array<char, Caller::KmerSize> &kmer, size_t offset) {
	// std::cerr << "process " << kmer << "\n";
	
	int ec;
	caller_.fetch(kmer, [this, hit_cb, offset, &hits, &calls, &otu_stats, &current_fI, &kmer, seqlen]
		      (const typename Caller::KData *kdata) {
	    /*
	      double z = kdata->var == 0 ? 0.0 : (seqlen - static_cast<double>(kdata->mean)) /
	      sqrt(static_cast<double>(kdata->var));
	      std::cerr << "fetch got " << kdata->function_index << " " << kdata->mean << " " << kdata->median << " " << kdata->var << " " << z << "\n";
	    */
	    if (hit_cb != nullptr)
		hit_cb(hit_in_sequence_t<Caller>({0L, kdata->otu_index, kdata->avg_from_end, kdata->function_index, kdata->weight}, offset, kmer, kdata));

	    // std::cerr << kmer << "\t" << offset << "\t" << kdata->function_index << "\n";
	    // Is this hit beyond max_gap_ of the last one?
	    if (!hits.empty() && hits.last_hit().pos + max_gap_ < offset)
	    {
		if (hits.count() >= min_hits_)
		    hits.process(current_fI, calls, otu_stats);
		else
		    hits.clear();
	    }
	    if (hits.empty())
	    {
		current_fI = kdata->function_index;
	    }
	    
	    if (!order_constraint_ ||
		hits.empty() ||
		(kdata->function_index == hits.last_hit().kdata.function_index &&
		 labs((offset - hits.last_hit().pos) -
		      (hits.last_hit().kdata.avg_from_end - kdata->avg_from_end)) <= 20))
	    {
//		using hittype = typename HitSet<Caller,K>::hit;
//		hits.emplace_back(hittype {*kdata, offset} );
		hits.emplace_back(typename HitSet<Caller>::hit{*kdata, offset} );
		/*
		 * If we have a pair of new functions, it is time to
		 * process one set and initialize the next.
		 */
		if (hits.count() > 1 && current_fI != kdata->function_index)
		{
		    auto end = hits.rbegin();
		    if (end[1].kdata.function_index == end[0].kdata.function_index)
		    {
			hits.process(current_fI, calls, otu_stats);
		    }
		}
	    }

	}, ec);
	if (ec)
	{
	    // std::cerr << "Error " << ec << " returned from caller\n";
	}
    });
    if (hits.count() >= min_hits_)
	hits.process(current_fI, calls, otu_stats);
}

template <class Caller>
void KmerGeneric<Caller>::process_aa_seq(const std::string &idstr, const std::string &seqstr,
				 std::shared_ptr<std::vector<KmerCall>> calls,
				 std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
				 std::shared_ptr<KmerOtuStats> otu_stats)
{
    gather_hits(seqstr, calls, hit_cb, otu_stats);

    if (otu_stats)
	otu_stats->finalize();
}

template <class Caller>
void KmerGeneric<Caller>::process_seq(const char *id,const char *data,
			   std::shared_ptr<std::vector<KmerCall>> calls,
			   std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
			   std::shared_ptr<KmerOtuStats> otu_stats)

{
}

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  

template <class Caller>
std::string KmerGeneric<Caller>::format_hit(const hit_in_sequence_t<Caller> &h)
{
    std::ostringstream oss;

    oss << "HIT\t" << h.offset << "\t" << "\t" << h.hit.avg_from_end << "\t" << function_at_index(h.hit.function_index) << "\t" << h.hit.weight << "\t" << h.hit.otu_index << "\n";
    
    return oss.str();
}


template <class Caller>
std::string KmerGeneric<Caller>::format_call(const KmerCall &c)
{
    std::ostringstream oss;
    oss << "CALL\t" << c.start << "\t" << c.end << "\t" << c.count;
    oss << "\t" << c.function_index << "\t" << function_at_index(c.function_index);
    oss << "\t" << c.weighted_hits
	//R2 << "\t" << c.r2
	<< "\n";

    return oss.str();
}


struct FuncScore
{
    int count;
    float weighted;
    FuncScore() : count(0), weighted(0.0) {}
    FuncScore(int c, float w) : count(c), weighted(w) {}
    FuncScore(const FuncScore &f) : count(f.count), weighted(f.weighted) {}
    FuncScore& operator=(const FuncScore& other)
	{
	    // check for self-assignment
	    if(&other == this)
		return *this;
	    count = other.count;
	    weighted = other.weighted;
	    return *this;
	}
};

/*
 * Find the best call from this set of calls.
 *
 * This code replicates the amino acid version of the SEED pipeline
 * km_process_hits_to_regions | km_pick_best_hit_in_peg
 */
template <class Caller>
void KmerGeneric<Caller>::find_best_call(const std::string &id, std::vector<KmerCall> &calls,
					 FunctionIndex &function_index, std::string &function,
					 float &score, float &weighted_score, float &score_offset, bool call_ambigs)
{
    function_index = UndefinedFunction;
    function = "";
    score = 0.0;
    weighted_score = 0.0;

    if (calls.size() == 0)
    {
	return;
    }
    
    /*
     * First merge adjacent hits that have the same function.
     */
    std::vector<KmerCall> collapsed;

    auto comp = calls.begin();

    while (comp != calls.end())
    {
	collapsed.push_back(*comp);
	comp++;
	KmerCall &cur = collapsed.back();

	while (comp != calls.end() && cur.function_index == comp->function_index)
	{
	    cur.end = comp->end;
	    cur.count += comp->count;
	    cur.weighted_hits += comp->weighted_hits;
	    comp++;
	}
    }
#if DEBUG_SCORING
    std::cout << "after collapse:\n";
    for (auto iter = collapsed.begin(); iter != collapsed.end(); iter++)
    {
	std::cout << format_call(*iter);
    }
#endif

    /*
     *
     * Merge hits when we have a case with
     * +------+--+-------+--+-----------+
     * |  F1  |  |   F2  |  |   F1      |
     * +------+--+-------+--+-----------+
     *
     * where the score for F2 is below 5 and the combined scores for
     * the two F1 hits is 10 or more.
     *
     * If that is the case we discard the F2 hit and combine
     * the F1 hits.
     */

    std::vector<KmerCall> merged;

    int merge_interior_thresh = 5;
    int merge_exterior_thresh = 10;

    comp = collapsed.begin();
    while (comp != collapsed.end())
    {
	merged.push_back(*comp);
	comp++;
	auto comp2 = comp + 1;
	KmerCall &cur = merged.back();
	while (comp != collapsed.end() && comp2 != collapsed.end() &&
	       cur.function_index == comp2->function_index &&
	       comp->count < merge_interior_thresh &&
	       (cur.count + comp2->count) >= merge_exterior_thresh)
	{
	    cur.end = comp2->end;
	    cur.count += comp2->count;
	    cur.weighted_hits += comp2->weighted_hits;
	    comp += 2;
	    comp2 = comp + 1;
	}
    }
   
#if DEBUG_SCORING
    std::cerr << "after merge:\n";
    for (auto iter = merged.begin(); iter != merged.end(); iter++)
    {
	std::cerr << format_call(*iter);
    }
#endif

    /*
     * Determine best call.
     *
     * In the current perl kmer_search (km_pick_best_hit_in_peg) code we just take the best
     * function in terms of weighted score. However, that allows tied scores to be called
     * arbitrarily, and close-to-tied scores to be settled based on insignificant differences.
     *
     * In the original kmerv1 code, we required a score threshold difference (typically 5 hits)
     * between the best function and next best function. We resurrect that here.
     *
     */

    /*
     * Fusions
     *
     * We make an attempt to call fusions here.
     */

    if (merged.size() > 1)
    {
	char next_func_key = 'A';
	char next_fusion_key = 'W';
	std::map<std::string, char> func_map;
	std::map<std::string, char> fusion_map;
	std::map<char, std::pair<FunctionIndex, std::string>> key_to_function_info;
	namespace acc = boost::accumulators;
	using accum = acc::accumulator_set<float, acc::stats<acc::tag::mean, acc::tag::variance>>;
	std::map<char, accum> part_stats;
	std::string exp;

	struct fusion_part
	{
	    char ident;
	    const KmerCall call;
	};
	std::vector<fusion_part> fusion_parts;
	float sum_weights = 0.0;
	int sum_scores = 0;
	for (auto c: merged)
	{
	    /* accumulate scores in the event we call the fusion */
	    sum_weights += c.weighted_hits;
	    sum_scores += c.count;
	       
	    auto func = function_at_index(c.function_index);
	    std::vector<std::string> parts = split(func, " / ");
	    std::string fusion_key;
	    for (auto part: parts)
	    {
		if (func_map.find(part) == func_map.end())
		{
		    char f = next_func_key++;
		    func_map[part] = f;
		}
		fusion_key += func_map[part];
	    }
	    if (parts.size() > 1)
	    {
		if (fusion_map.find(fusion_key) == fusion_map.end())
		{
		    fusion_map[fusion_key] = next_fusion_key++;
		}
		char fkey = fusion_map[fusion_key];
		exp += fkey;
		fusion_parts.emplace_back(fusion_part { fkey, c });
		part_stats[fkey](static_cast<float>(c.protein_length_median));

		key_to_function_info[fkey] = std::make_pair(c.function_index, func);
	    }
	    else
	    {
		exp += func_map[func];
		fusion_parts.emplace_back(fusion_part { func_map[func], c });
		part_stats[func_map[func]](static_cast<float>(c.protein_length_median));
		key_to_function_info[func_map[func]] = std::make_pair(c.function_index, func);
	    }
		
	}
//	std::cerr << id << "\t" << exp << "\n";
#if DEBUG_SCORING
	std::cerr << "Exp list: " << exp << "\n";
	for (auto x: fusion_parts) {
	    std::cerr << x.ident << ": " << x.call << "\n";
	}
#endif

	const boost::regex fusion_re("^W?A[A|W]*W[B|W]*BW?");
	if (regex_match(exp, fusion_re))
	{
	    #if DEBUG_SCORING
	    std::cerr << "potential fusion " << exp << "\n";
	    for (auto p: { 'A', 'W', 'B' })
	    {
		auto stats = part_stats[p];
		std::cerr << p << " " << acc::mean(stats) << " " << acc::variance(stats) << "\n";
	    }
	    #endif
	    auto a_mean = acc::mean(part_stats['A']);
	    auto w_mean = acc::mean(part_stats['W']);
	    auto b_mean = acc::mean(part_stats['B']);
	    auto diff = (a_mean + b_mean) - w_mean;
	    auto frac_dif = abs(diff) / w_mean;
	    if (frac_dif < 0.1)
	    {
		/*
		  std::cerr << "call fusion " << id << " " << exp << "\n";
		std::cerr << a_mean << " " << w_mean << " " << b_mean << " " << diff  << " " << frac_dif << "\n";
		for (auto x: fusion_parts) 
		    std::cerr << x.ident << ": " << x.call << "\n";
		*/
		function_index = key_to_function_info['W'].first;
		function = key_to_function_info['W'].second;
		score = sum_scores;
		weighted_score = sum_weights;
		score_offset = 0.0;
		return;
	    }
	    else
	    {
		/*
		std::cerr << "no-call fusion " << id << " " << exp << "\n";
		std::cerr << a_mean << " " << w_mean << " " << b_mean << " " << diff  << " " << frac_dif << "\n";
		*/
	    }
	}
	    
    }

    typedef std::map<int, FuncScore> map_t;

    map_t by_func;

    for (auto c: merged)
    {
	auto it = by_func.find(c.function_index);
	if (it == by_func.end())
	{
	    by_func.insert(std::make_pair(c.function_index, FuncScore(c.count, c.weighted_hits)));
	}
	else
	{
	    it->second.count += c.count;
	    it->second.weighted += c.weighted_hits;
	}
    }

    //typedef map_t::value_type ent_t;
    typedef std::pair<FunctionIndex, FuncScore> ent_t;

    std::vector<ent_t> vec;
    for (auto it = by_func.begin(); it != by_func.end(); it++)
	vec.push_back(*it);

//    std::cerr << "vec len " << vec.size() << "\n";
    if (vec.size() > 1)
    {
	std::partial_sort(vec.begin(), vec.begin() +  2, vec.end(), 
			  [](const ent_t& s1, const ent_t& s2) {
			      return (s1.second.count > s2.second.count); });
//			      return (s1.second.weighted > s2.second.weighted) || (s1.second.count > s2.second.count); });
    }
    
#if DEBUG_SCORING
    for (auto x: vec)
    {
	std::cerr << x.first << " " << x.second.count << " " << x.second.weighted << " ";
	std::cerr << function_at_index(x.first) << "\n";
    }
#endif
    
    if (vec.size() == 1)
	score_offset = (float) vec[0].second.count;
    else
	score_offset = (float) (vec[0].second.count - vec[1].second.count);

#if DEBUG_SCORING
    std::cerr << "Offset=" << score_offset << "\n";
#endif

    if (score_offset >= 5.0)
    {
	auto best = vec[0];
	function_index = best.first;
	function = function_at_index(function_index);
	score = (float) best.second.count;
	weighted_score = best.second.weighted;
    }
    else
    {
	function_index = UndefinedFunction;
	function = "";
	score = 0.0;

	/*
	 * Try to compute a fallback function naming the two best hits if there are two hits within the
	 * threshold but greater than the next hit.
	 */
	if (call_ambigs)
	{
	    if (vec.size() >= 2)
	    {
		std::string f1 = function_at_index(vec[0].first);
		std::string f2 = function_at_index(vec[1].first);
		if (f2 > f1)
		    std::swap(f1, f2);
		
		if (vec.size() == 2)
		{
		    function = f1 + " ?? " + f2;
		    score = (float) vec[0].second.count;
		}
		else if (vec.size() > 2)
		{
		    float pair_offset = (float) (vec[1].second.count - vec[2].second.count);
		    if (pair_offset > 2.0)
		    {
			function = f1 + " ?? " + f2;
			score = (float) vec[0].second.count;
			score_offset = pair_offset;
			weighted_score = vec[0].second.weighted;
		    }
		}
	    }
	}		
    }
}

template <class Caller>
template <typename HitCB, typename CallCB>
void KmerGeneric<Caller>::process_fasta_stream(std::istream &istr, HitCB &hit_cb, CallCB &call_cb)
{

    try {
	FastaParser parser;
    
	parser.set_callback([this, &hit_cb, &call_cb](const std::string &id, const std::string &seq) {
	
	    double slen = static_cast<double>(seq.length());
	    auto calls = std::make_shared<std::vector<KmerCall>>();
	    auto ostats = std::make_shared<KmerOtuStats>();

	    /*
	    auto hit_cb = [ slen, &this](const hit_in_sequence_t<Caller> &hit) {
		auto kd = hit.kdata;
		if (false)
		{
		    std::cerr << hit.kmer << "\t" << function_at_index(kd->function_index) << "\t" << kd->median << "\t" << kd->mean << "\t" << kd->var << "\t" << sqrt(kd->var) << "\t" << "\n";
		}
		};
	    */

	    // std::cerr << "stream " << id << " " << seq << "\n";
	    process_aa_seq(id, seq, calls, hit_cb, ostats);
	    for (auto c: *calls)
	    {
		// std::cerr << c << "\n";
	    }
	    FunctionIndex fi;
	    std::string func;
	    float score;
	    float wt;
	    float offset;
	    find_best_call(id, *calls, fi, func, score, wt, offset);
	    call_cb(id, func, fi, score);
	    // std::cerr << id << "\t" << func << "\t" << fi << "\t" << score << "\n";
	    return 0;
	});

	parser.parse(istr);
	parser.parse_complete();
    }
    catch (std::runtime_error &x)
    {
	std::cerr<< "caught " << x.what() << "\n";
    }
}
