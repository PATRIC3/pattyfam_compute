#ifndef _fasta_parser_h
#define _fasta_parser_h

#include <sstream>
#include <iostream>
#include <functional>

class FastaParser
{
public:

    enum state {
	s_start = 0,
	s_id,
	s_defline,
	s_data,
	s_id_or_data,
    };

    FastaParser();

    void set_callback(std::function<void(const std::string &id, const std::string &seq)> on_seq) {
	on_seq_ = on_seq;
    }
    
    void set_def_callback(std::function<void(const std::string &id, const std::string &def, const std::string &seq)> on_seq) {
	on_def_seq_ = on_seq;
    }

    void set_error_callback(std::function<bool(const std::string &err, int line, const std::string id)> cb)
    {
	on_error_ = cb;
    }

    void parse(std::istream &stream);
    void init_parse();
    // return true to continue parsing
    inline bool parse_char(char c)
    {
	/*
	  std::cout << "top " << c << " " << (0+cur_state_) << std::endl;
	  std::cout << cur_id_ << "\n";
	  std::cout << cur_seq_ << "\n";
	*/
	if (c == '\n')
	    line_number++;
	if (c == '\r')
	    return true;
	
	std::string err;
	switch (cur_state_)
	{
	case s_start:
	    if (c != '>')
	    {
		err = "Missing >";
	    }
	    else
	    {
		cur_state_ = s_id;
	    }
	    break;
	    
	case s_id:
	    if (isblank(c))
	    {
		cur_def_.push_back(c);
		cur_state_ = s_defline;
	    }
	    else if (c == '\n')
	    {
		cur_state_ = s_data;
	    }
	    else
	    {
		cur_id_.push_back(c);
	    }
	    break;

	case s_defline:
	    if (c == '\n')
	    {
		cur_state_ = s_data;
	    }
	    else
	    {
		cur_def_.push_back(c);
	    }
	    break;

	case s_data:
	    if (c == '\n')
	    {
		cur_state_ = s_id_or_data;

	    }
	    else if (isalpha(c) || c == '*')
	    {
		cur_seq_.push_back(c);
	    }
	    else
	    {
		err = "Bad data character '";
		err += c;
		err += "'"; 
	    }
	    break;

	case s_id_or_data:
	    if (c == '>')
	    {
		call_callback();
		cur_id_ = "";
		cur_def_ = "";
		cur_seq_ = "";
		cur_state_ = s_id;
	    }
	    else if (c == '\n')
	    {
		// No change. cur_state_ = s_id_or_data;
	    }
	    else if (isalpha(c))
	    {
		cur_seq_.push_back(c);
		cur_state_ = s_data;
	    }
	    else
	    {
		err = "Bad id or data character '";
		err += c;
		err += "'";
	    }
	    break;
	}
	if (!err.empty())
	{
	    std::cerr << "Error found: " << err << " at line " << line_number << " id='" << cur_id_ << "'" << std::endl;
	    if (on_error_)
	    {
		return on_error_(err, line_number, cur_id_);
	    }
	}
	return true;
    }

    void parse_complete();

private:
    int line_number;
    state cur_state_;
    std::string cur_id_;
    std::string cur_seq_;
    std::string cur_def_;
    
    std::function<void(const std::string &id, const std::string &seq)> on_seq_;
    std::function<void(const std::string &id, const std::string &def, const std::string &seq)> on_def_seq_;
    std::function<bool(const std::string &err, int line, const std::string id)> on_error_;
    void call_callback()
    {
	if (on_seq_)
	    on_seq_(cur_id_, cur_seq_);
	if (on_def_seq_)
	    on_def_seq_(cur_id_, cur_def_, cur_seq_);
    }
};


#endif
