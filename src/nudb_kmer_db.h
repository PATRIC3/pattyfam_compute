#ifndef _nudb_kmer_db_h
#define _nudb_kmer_db_h

/**
 * Kmer database using NuDB.
 *
 * Encapsulate all the database ops here.
 */

#include <nudb/nudb.hpp>
#include <iostream>
#include <filesystem>
#include <string>
#include <array>
#include <cstdint>

#include "kmer_value_types.h"

template <int K>
class NuDBKmerDb
{
public:
    static const int KmerSize = K;
    struct KData {
	OTUIndex otu_index = UndefinedOTU;
	uint16_t  avg_from_end = 0;
	FunctionIndex function_index = UndefinedFunction;
	float weight = 0.0;

#if 0
	unsigned int mean = 0;
	unsigned int median = 0;
	unsigned int var = 0;
#else
	
	unsigned short mean = 0;
	unsigned short median = 0;
	unsigned short var = 0;
	unsigned short n_proteins = 0;
#endif
    };
    using key_type = std::array<char, K>;
    static constexpr int kmer_size = K;

    NuDBKmerDb(const std::string &file_base)
	: file_base_(file_base)
	, dat_path_(file_base + ".dat")
	, key_path_(file_base + ".key")
	, log_path_(file_base + ".log") {
    }

    ~NuDBKmerDb() {

	if (db_.is_open())
	{
	    nudb::error_code ec;
	    // std::cerr << "close db\n";
	    db_.close(ec);
	    if (ec)
		std::cerr << "Error during db close: " << ec.message() << "\n";
	}
    }

    nudb::store &db() { return db_; }

    bool exists() {
	return std::filesystem::exists(dat_path_);
    }

    void create(nudb::error_code &ec){
	nudb::create<nudb::xxhasher>(dat_path_, key_path_, log_path_,
				     1,
				     nudb::make_salt(),
				     sizeof(key_type),
				     nudb::block_size("."),
				     0.5f,
				     ec);
    }

    void create() {
	nudb::error_code ec;
	create(ec);
	if (ec)
	{
	    throw std::runtime_error("error in NuDBKmerDb::create(): " + ec.message());
	}
    }

    void open(nudb::error_code &ec) {
	db_.open(dat_path_, key_path_, log_path_, ec);
    }

    void open() {
	nudb::error_code ec;
	open(ec);
	if (ec)
	{
	    throw std::runtime_error("error in NuDBKmerDb::open(): " + ec.message());
	}
    }

    void insert(const std::string &key, const KData &kdata, nudb::error_code &ec) {
	key_type ka;
	if (key.length() != kmer_size)
	    throw std::runtime_error("Invalid kmer size");
	std::copy(key.begin(), key.end(), ka.data());
	insert(ka, kdata, ec);
    }
    void insert(const key_type &key, const KData &kdata, nudb::error_code &ec) {
	db_.insert(key.data(), &kdata, sizeof(kdata), ec);
    }

    template <typename CB>
    void fetch(const key_type &key, CB cb, nudb::error_code &ec) {
	db_.fetch(key.data(), [&cb](void const *buffer,  std::size_t size) {
	    if (size != sizeof(KData))
	    {
		std::cerr << "Invalid data size: " << size << " != " << sizeof(KData) << "\n";
	    }
	    assert(size == sizeof(KData));
	    const KData *kdata = static_cast<const KData *>(buffer);
	    cb(kdata);
	}, ec);
    }
    template <typename CB>
    void fetch(const key_type &key, CB cb, int &iec) {
	nudb::error_code ec;
	db_.fetch(key.data(), [&cb](void const *buffer,  std::size_t size) {
	    if (size != sizeof(KData))
	    {
		std::cerr << "Invalid data size: " << size << " != " << sizeof(KData) << "\n";
	    }
	    assert(size == sizeof(KData));
	    const KData *kdata = static_cast<const KData *>(buffer);
	    cb(kdata);
	}, ec);
	// std::cerr << ec.message() << std::endl;
	iec = ec.value();
    }

private:
    std::string file_base_;
    std::string dat_path_, key_path_, log_path_;
    nudb::store db_;
    
};


#endif // _nudb_kmer_db_h
