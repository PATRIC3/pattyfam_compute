#ifndef _OPERATORS_H
#define _OPERATORS_H

#include <functional>
#include <iostream>
#include <string>

template<class T>
struct less_second : std::binary_function<T,T,bool>
{
    inline bool operator()( const T& lhs, const T& rhs )
	{
	    return rhs.second < lhs.second;
	}
};  

template<typename T>
class PrefexOutputIterator
{
    std::ostream&       ostream;
    std::string         prefix;
    bool                first;
public:

    typedef std::size_t                 difference_type;
    typedef T                           value_type;
    typedef T*                          pointer;
    typedef T                           reference;
    typedef std::output_iterator_tag    iterator_category;

PrefexOutputIterator(std::ostream& o,std::string const& p = ""): ostream(o), prefix(p), first(true) {}

    PrefexOutputIterator& operator*()       {return *this;}
    PrefexOutputIterator& operator++()      {return *this;}
    PrefexOutputIterator& operator++(int)   {return *this;}

    void operator=(T const& value)
    {
	if (first)      {ostream << value;first = false;}
	else            {ostream << prefix << value;}
    }
};

template <class MapType>
class MapKeyIterator {
public:
    class iterator {
    public:
    iterator(typename MapType::iterator it) : it(it) {}
	iterator operator++() { return ++it; }
	bool operator!=(const iterator & other) { return it != other.it; }
	typename MapType::key_type operator*() const { return it->first; }  // Return key part of map
    private:
	typename MapType::iterator it;
    };
private:
    MapType& map;
public:
MapKeyIterator(MapType& m) : map(m) {}
    iterator begin() { return iterator(map.begin()); }
    iterator end() { return iterator(map.end()); }
};
template <class MapType>
MapKeyIterator<MapType> MapKeys(MapType& m)
{
    return MapKeyIterator<MapType>(m);
}


template <typename K, typename V>
    std::vector<std::pair<K, V>> sort_by_values(std::map<K, V> &map)
{
    std::vector<std::pair<K, V>> vec;
    for (auto x: map)
	vec.push_back(x);
    std::sort(vec.begin(), vec.end(), less_second<std::pair<K, V>>());
    return vec;
}

std::vector<std::string> split(const std::string &s, const std::string &delim){
    std::vector<std::string> result;
    int start = 0;
    std::string::size_type end = 0;
    while(end != std::string::npos)
    {
	end = s.find(delim, start);
	result.push_back(s.substr(start, end-start));
	start = end + delim.length();
    }
    return result;
}

template <typename Elt>
std::ostream &operator<<(std::ostream &os, std::vector<Elt> &v)
{
    bool first = true;
    for (auto x: v)
    {
	if (!first)
	{
	    os << " ";
	} else {
	    first = false;
	}
	os << x;
    }
    return os;
}
		

#endif
