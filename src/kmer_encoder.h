#ifndef _kmer_encoder_h
#define _kmer_encoder_h

#include <cstdint>
#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>

#include "kmer_params.h"

class KmerEncoder
{
public:
    KmerEncoder();

    static struct AAText_t  {} AAText;
    static struct AANumeric_t {} AANumeric;

    inline unsigned char to_amino_acid_off(uint8_t c) {
	return aa_to_offset_[c];
    };

    unsigned long long encoded_kmer(const std::string &p) {
	return encoded_kmer(reinterpret_cast<const unsigned char *>(p.c_str()));
    }

    unsigned long long encoded_kmer(const unsigned char *p) {
	// std::cout << "encode '" << p << "'\n";
	unsigned long long encodedK = *p;
	int i;

	for (i=1; (i <= KMER_SIZE-1); i++) {
	    encodedK = (encodedK * 20) + p[i];
	}

	if (encodedK > MAX_ENCODED) {
	    fprintf(stderr,"bad encoding - input must have included invalid characters\n");
	    for (i=0; (i < KMER_SIZE); i++) {
		fprintf(stderr,"%d ",*(p+i));
	    }
	    fprintf(stderr,"\n");
	    exit(2);
	}
	return encodedK;
    }

    unsigned long long encoded_aa_kmer(const std::string &p) {
	return encoded_aa_kmer(p.c_str());
    }

    unsigned long long encoded_aa_kmer(const unsigned char *p) {
	return encoded_aa_kmer(reinterpret_cast<const char *>(p));
    }

    unsigned long long encoded_aa_kmer(const char *p) {
	unsigned char aa_off[KMER_SIZE];
	int j;
	for (j=0; (j < KMER_SIZE); j++) {
	    unsigned char prot_c = (unsigned char) p[j];
	    aa_off[j] = to_amino_acid_off(prot_c);
	    /*
	     * If we have an invalid char, return an invalid encoding.
	     */
	    if (aa_off[j] >= 20)
		return MAX_ENCODED + 1;
	}
	return encoded_kmer(aa_off);
    }

    template <class T, size_t N>
	unsigned long long encoded_aa_kmer(std::array<T, N> &p) {
	unsigned char aa_off[N+1];
	unsigned char *ap = aa_off;
	for (auto prot_c: p)
	{
	    *ap = to_amino_acid_off(prot_c);
	    /*
	     * If we have an invalid char, return an invalid encoding.
	     */
	    if (*ap >= 20)
		return MAX_ENCODED + 1;
	    ap++;
	}
	*ap = 0;
	return encoded_kmer(aa_off);
    }

    void decoded_kmer(unsigned long long encodedK,char *decoded) {
  
	int i;
	*(decoded+KMER_SIZE) = '\0';
	unsigned long long x = encodedK;
	
	for (i=KMER_SIZE-1; (i >= 0); i--) {
	    *(decoded+i) = prot_alpha[x % 20];
	    x = x / 20;
	}
    }

    /**
     * Translate a dna sequence (starting at offset off in the string)
     * to protein; we also create the encoded numeric amino acid number.
     */
    void translate(const std::string &seq, int off, std::string &aa, std::string &aa_num) {
	aa.clear();
	aa.reserve(seq.size() / 3 + 2);
	aa_num.clear();
	aa_num.reserve(seq.size() / 3 + 2);

	size_t max = seq.size() - 3;

	for (size_t i = (size_t) off; (i <= max); )
	{
	    while (std::isspace(seq[i]))
		++i;
	    int c1 = dna_char(seq[i++]);
	    while (std::isspace(seq[i]))
		++i;
	    int c2 = dna_char(seq[i++]);
	    while (std::isspace(seq[i]))
		++i;
	    int c3 = dna_char(seq[i++]);
	    if ((c1 < 4) && (c2 < 4) && (c3 < 4)) {
		int I = (c1 * 16) + (c2 * 4) + c3;
		char prot_c = genetic_code[I];
		aa += prot_c;
		aa_num += to_amino_acid_off(prot_c);
	    }
	    else {
		aa += 'x';
		aa_num += 20;
	    }
	}
	aa_num += 21;
    }

    void rev_comp(const std::string &fwd, std::string &rev) {
	for (auto it = fwd.rbegin(); it != fwd.rend(); it++)
	{
	    if (!isspace(*it))
		rev += complement_base(*it);
	}
    }

    template <int N, typename F>
    void for_each_kmer(const std::string &str, F cb) {
	const char *ptr = str.c_str();
	const char *end = ptr  + str.length();
	const char *last_kmer= end - N + 1;

	const char *next_ambig = std::find_if(ptr, end, [](char c) -> bool { return c == '*' || c == 20; });
	
	char buf[N+1];
	while (ptr <= last_kmer)
	{
	    const char *kend = ptr + N - 1;
	    if (kend >= next_ambig)
	    {
		ptr = next_ambig + 1;
		next_ambig = std::find_if(ptr, end, [](char c) { return c == '*' || c == 20; });
		continue;
	    }
	    strncpy(buf, ptr, N);
	    buf[N] = 0;
	    cb(buf);
	    ptr++;
	}
    }
    template <int N, typename F>
    void for_each_encoded_kmer(AANumeric_t, const std::string &str, F cb) {
	const unsigned char *start = reinterpret_cast<const unsigned char *>(str.c_str());
	const unsigned char *ptr = start;
	const unsigned char *end = ptr  + str.length();
	const unsigned char *last_kmer= end - N + 1;

	auto is_ambig = [](char c) -> bool { return c == 20; };

	const unsigned char *next_ambig = std::find_if(ptr, end, is_ambig);
	while (ptr <= last_kmer)
	{
	    const unsigned char *kend = ptr + N - 1;
	    if (kend >= next_ambig)
	    {
		ptr = next_ambig + 1;
		next_ambig = std::find_if(ptr, end, is_ambig);
		continue;
	    }
	    auto enc = encoded_kmer(ptr);
	    cb(ptr - start, enc);
	    ptr++;
	}
    }

    template <int N, typename F>
    void for_each_encoded_kmer(AAText_t, const std::string &str, F cb) {
	const unsigned char *start = reinterpret_cast<const unsigned char *>(str.c_str());
	const unsigned char *ptr = start;
	const unsigned char *end = ptr  + str.length();
	const unsigned char *last_kmer= end - N + 1;

	auto is_ambig = [](char c) -> bool { return c == '*' || c == 'X'; };

	const unsigned char *next_ambig = std::find_if(ptr, end, is_ambig);
	while (ptr <= last_kmer)
	{
	    const unsigned char *kend = ptr + N - 1;
	    if (kend >= next_ambig)
	    {
		ptr = next_ambig + 1;
		next_ambig = std::find_if(ptr, end, is_ambig);
		continue;
	    }
	    auto enc = encoded_aa_kmer(ptr);
	    cb(ptr - start, enc);
	    ptr++;
	}
    }

    char complement_base(char c)
	{
	    switch (c)
	    {
	    case 'a':
		return 't';
	    case 'A':
		return 'T';

	    case 'c':
		return 'g';
	    case 'C':
		return 'G';

	    case 'g':
		return 'c';
	    case 'G':
		return 'C';

	    case 't':
	    case 'u':
		return 'a';
	    case 'T':
	    case 'U':
		return 'A';

	    case 'm':
		return 'k';
	    case 'M':
		return 'K';

	    case 'r':
		return 'y';
	    case 'R':
		return 'Y';

	    case 'w':
		return 'w';
	    case 'W':
		return 'W';

	    case 's':
		return 'S';
	    case 'S':
		return 'S';

	    case 'y':
		return 'r';
	    case 'Y':
		return 'R';

	    case 'k':
		return 'm';
	    case 'K':
		return 'M';

	    case 'b':
		return 'v';
	    case 'B':
		return 'V';

	    case 'd':
		return 'h';
	    case 'D':
		return 'H';

	    case 'h':
		return 'd';
	    case 'H':
		return 'D';

	    case 'v':
		return 'b';
	    case 'V':
		return 'B';

	    case 'n':
		return 'n';
	    case 'N':
		return 'N';

	    default:
		return c;
	    }
	}


    int dna_char(char c)
	{
	    switch (c)
	    {
	    case 'a':
	    case 'A':
		return 0;

	    case 'c':
	    case 'C':
		return 1;

	    case 'g':
	    case 'G':
		return 2;

	    case 't':
	    case 'u':
	    case 'T':
	    case 'U':
		return 3;

	    default:
		return 4;;
	    }
	}

    uint8_t aa_to_offset_[256];
    static constexpr char prot_alpha[20] = {
	'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' 
    };

    static constexpr  char genetic_code[64] = {
    'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
    'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
    'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
    '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'
    };
};


#endif // _kmer_encoder_h
