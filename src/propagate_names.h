#ifndef _PROPAGATE_NAMES_H
#define _PROPAGATE_NAMES_H

#include <string>
#include <set>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <thread>
#include <boost/filesystem.hpp>
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_queue.h"

typedef tbb::concurrent_unordered_set<std::string> string_set_t;
typedef tbb::concurrent_unordered_map<std::string, std::string> string_map_t;
typedef tbb::concurrent_unordered_map<std::string, std::vector<std::string>> string_vec_map_t;
typedef tbb::concurrent_unordered_map<std::string, std::set<std::string>> string_set_map_t;
typedef tbb::concurrent_unordered_multimap<std::string, std::string> string_multimap_t;

class FamData
{
public:
    enum FamilyType
    {
	FAM_LOCAL,
	FAM_GLOBAL
    };

    FamData(const boost::filesystem::path &fams_file, const boost::filesystem::path &data_dir, const std::string &target_genus,
	    FamilyType family_type) :
    fams_file_(fams_file), data_dir_(data_dir), target_genus_(target_genus), family_type_(family_type)
    {
	int sz =  1000000;
	int sz2 =   1000000;
	fid_is_key_.rehash(sz);
	md5_to_key_.rehash(sz);
	fid_to_md5_.rehash(sz);
	fam_to_md5s_.rehash(sz);
	fam_to_function_.rehash(sz);
	md5_to_fam_.rehash(sz);
    }

    ~FamData()
    {
	std::cout << "Destroy famdata for " << data_dir_ << "\n";
    }
    
    void read_fams_file();

    void read_pegsyn_file(boost::filesystem::path &data_file);
    void read_pegsyn();

    bool exists(const std::string &md5)
    {
	return md5_to_key_.find(md5) != md5_to_key_.end();
    }

    const std::string peg_to_fam(const std::string &md5)
    {
	auto iter = md5_to_fam_.find(md5);
	if (iter == md5_to_fam_.end())
	{
	    std::cerr << "No family found for " << md5 << "\n";
	    return "";
	}
	else
	{
	    return iter->second;
	}
    }
    const std::string fam_to_fun(const std::string &fam)
    {
	auto iter = fam_to_function_.find(fam);
	if (iter == fam_to_function_.end())
	{
	    std::cerr << "No function found for " << fam << "\n";
	    return "";
	}
	else
	{
	    return iter->second;
	}
    }
		
    boost::filesystem::path fams_file_;
    boost::filesystem::path data_dir_;
    std::string target_genus_;
    FamilyType family_type_;

    // peg.synonyms support
    string_map_t fid_is_key_;
    string_map_t md5_to_key_;
    string_map_t fid_to_md5_;

    // family file support
    string_set_map_t fam_to_md5s_;
    string_map_t fam_to_function_;
    string_map_t md5_to_fam_;
};

std::ostream &operator<<(std::ostream &os, const string_set_t &d)
{
    for (auto x: d)
    {
	os << "   " << x << "\n";
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const string_map_t &d)
{
    for (auto x: d)
    {
	os << "   " << x.first << ": " << x.second << "\n";
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const string_multimap_t &d)
{
    for (auto x: d)
    {
	os << "   " << x.first << ": " << x.second << "\n";
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const string_vec_map_t &d)
{
    for (auto x: d)
    {
	os << "   " << x.first << ":";
	for (auto y: x.second)
	{
	    os << " " << y;
	}
	os << "\n";
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const string_set_map_t &d)
{
    for (auto x: d)
    {
	os << "   " << x.first << ":";
	for (auto y: x.second)
	{
	    os << " " << y;
	}
	os << "\n";
    }
    return os;
}


std::ostream &operator<<(std::ostream &os, const FamData &d)
{
    os << "fid_is_key:\n" << d.fid_is_key_ << "\n";
    os << "md5_to_key:\n" << d.md5_to_key_ << "\n";
    os << "fid_to_md5:\n" << d.fid_to_md5_ << "\n";
    os << "fam_to_md5s:\n" << d.fam_to_md5s_ << "\n";
    os << "fam_to_function:\n" << d.fam_to_function_ << "\n";
    os << "md5_to_fam:\n" << d.md5_to_fam_ << "\n";
    return os;
}

class RenumberState
{
public:
    RenumberState(FamData &old_data, FamData &new_data);

    FamData &old_data_;
    FamData &new_data_;

    void phase_1();
    void phase_1_body(const std::string &fam, std::set<std::string> &fids);

    void phase_2();
    void phase_2_body(const std::string &nfam, std::set<std::string> &nfids);

    void phase_3();
    void phase_3_body(const std::string &fam, std::set<std::string> &fids);

    void write_unmapped();

    void set_logfile(const std::string &file);

    std::string allocate_new_id();
    void log_result(const std::string &res);

    std::ofstream log_;

    tbb::concurrent_vector<std::string> results_;

    tbb::concurrent_unordered_map<std::string, std::set<std::string>> old_fam_to_new_fam_set_;

    string_map_t old_fam_used_;
    string_map_t new_fam_name_;

    tbb::atomic<int> new_idx_;

    void logger_main();
    void start_logger();
    void stop_logger();
    tbb::concurrent_bounded_queue<std::string> result_queue_;
    std::thread result_logger_thread_;
};

#endif
