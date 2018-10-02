#ifndef GENE
#define GENE
#include <set>
#include <genome.hpp>
#include <indexer.hpp>

struct inv {
  std::string chr;
  unsigned start, end;
  bool strand;
  std::string seq;
  bool operator==(const inv &another) {
    return ((chr == another.chr) && (start == another.start) &&
            (end == another.end) && (strand == another.strand));
  };
  bool operator!=(const inv &another) {
    return !((*this)==(another));
  };
};

struct interval {
	// note that intervals are "unsigned"
	// and only count absolute pos on + strand

public:
	interval(seqdb * sdb, const unsigned& chr_, const unsigned & start, const unsigned & end,
		const unsigned & idx, bool strand) : sdb(sdb), chr_(chr_), start(start), end(end),
		idx(idx), strand(strand) {};
	interval() {};
	unsigned start = 0, end = 0, idx = 0, chr_=0;
	bool strand;
    std::string chr() const {return sdb->chrs[chr_];};
    std::string ref() const {return sdb->name;};
    unsigned long sz() const {return sdb->sizes[chr()];};
	std::string seq(int offset = 0)
	{
		if (offset)
			return sdb->get(chr(), start, start + offset);
		else
			return sdb->get(chr(), start, end);
	};
	inv v() { return inv{ chr(), start, end, strand, seq() }; };
	bool empty() { return (start | end | idx); };
	std::string info()
	{
		char tmp[256];
		auto fm = sprintf(tmp,
			"[%s, %s]: (%d, %d) @ %s @ {index: %d , strand: %d}",
			ref().c_str(), chr().c_str(), start, end, seq(50).c_str(), idx, (int)strand);
		return std::string(tmp);
	};
	bool operator==(const interval &another)
	{
		return ((chr() == another.chr()) && (start == another.start) &&
			(end == another.end) && (strand == another.strand));
	};
	bool operator!=(const interval &another)
	{
		return !((*this) == (another));
	};
private:
	seqdb *sdb = nullptr;
};

struct intervalCmp {
	bool operator()(const interval &l, const interval &r) {
		if (l.chr() != r.chr())
			return l.chr() < r.chr();
		if (r.start > l.end)
			return true;
		if (r.end < l.start)
			return false;
		return (l.start < r.start);
	};
};

bool combine_intervals(interval &left, interval &right) {
  if (left.ref().find(right.ref()) == std::string::npos)
    return false;
  if (left.chr().find(right.chr()) == std::string::npos)
    return false;
  if (left.start > right.start)
    left.start = right.start;
  if (left.end < right.end)
    left.end = right.end;
  return true;
}
class gene {
  private:
    seqdb* sdb = nullptr;
  public:
    std::string sym, id = "";
    unsigned cds_start, cds_end, tx_start, tx_end, idx, chr_;
    bool strand;
    std::multiset<interval, intervalCmp> exons;
    interval inv() { return interval{sdb, chr_, tx_start, tx_end, idx, strand}; };
    interval cds() { return interval{sdb, chr_, cds_start, cds_end, idx, strand}; };
    std::string chr() const {return sdb->chrs[chr_];};
    std::string ref() const {return sdb->name;};
    unsigned long sz() const {return sdb->sizes[chr()];};
    std::string info()
    {
        return inv().info();
    };
    std::vector<interval> get_exons()
    {
        std::vector<interval> out;
        for (auto &it : exons)
        {
            out.push_back(it);
        }
        return out;
    };
    std::vector<interval> get_introns()
    {
        std::vector<interval> out;
        if (exons.size() < 2)
            return out; // no intron in this situation
        auto it = exons.begin();
        auto itt = ++it;
        --it;
        while (itt != exons.end())
        {
            out.push_back(interval{sdb, chr_, it->end+1, itt->start-1, idx, strand});
            ++it;
            ++itt;
        }
        return out;
    };
    interval utr(const bool s)
    {
        unsigned start, end;
        if (s)
        {
            start = tx_start;
            end = cds_start - 1;
        }
        else
        {
            start = cds_end + 1;
            end = tx_end;
        }
        return interval{sdb, chr_, start, end, idx, strand};
    };
    interval utr5()
    {
        if (noncoding())
        {
            //std::cerr << "Gene at " << idx << " is noncoding, no 5'UTR" << std::endl;
            return interval();
        }
        return utr(strand);
    };
    interval utr3()
    {
        if (noncoding())
        {
            //std::cerr << "Gene at " << idx << " is noncoding, no 3'UTR" << std::endl;
            return interval();
        }
        return utr(!strand);
    };
    interval get_p(const int &l, const int &r, const bool &s)
    {
        unsigned start, end;
        if (s)
        {
            start = tx_start - l;
            end = tx_start + r;
        }
        else
        {
            start = tx_end + l;
            end = tx_end - r;
        }
        return interval{sdb, chr_, start, end, idx, strand};
    };
    interval get_promoter(const int &l, const int &r)
    {
        if (noncoding())
        {
            std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                      << std::endl;
            return interval();
        }
        return get_p(l, r, strand);
    };
    interval get_tail(const int &l, const int &r)
    {
        if (noncoding())
        {
            std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                      << std::endl;
            return interval();
        }
        return get_p(l, r, !strand);
    };
    virtual bool noncoding()
    {
        if (cds_start == cds_end)
        {
            return true;
        }
        return false;
    };
    bool operator==(const gene &another) { return (idx == another.idx); };
    bool operator!=(const gene &another) { return (idx != another.idx); };
};
#endif