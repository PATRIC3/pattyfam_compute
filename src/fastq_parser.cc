#include "fastq_parser.h"
#include <cctype>

FastqParser::FastqParser() : line_number(1)
{
    init_parse();
}

void FastqParser::init_parse()
{
    cur_state_ = s_start;
    cur_id_ = "";
    cur_def_ = "";
    cur_seq_ = "";
}

void FastqParser::parse(std::istream &stream)
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

void FastqParser::parse_complete()
{
    call_callback();
    cur_id_ = "";
    cur_def_ = "";
    cur_seq_ = "";
}

	
