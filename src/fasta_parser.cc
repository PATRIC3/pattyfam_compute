#include "fasta_parser.h"
#include <cctype>

FastaParser::FastaParser() : line_number(1)
{
    init_parse();
}

void FastaParser::init_parse()
{
    cur_state_ = s_start;
    cur_id_ = "";
    cur_def_ = "";
    cur_seq_ = "";
}

void FastaParser::parse(std::istream &stream)
{
    char c;
    init_parse();
    while (stream.get(c))
    {
	bool ok = parse_char(c);
	if (!ok)
	    break;
    }
    parse_complete();
}

void FastaParser::parse_complete()
{
    if (cur_id_[0])
	call_callback();
    cur_id_ = "";
    cur_def_ = "";
    cur_seq_ = "";
}

	
