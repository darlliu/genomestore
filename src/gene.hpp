#ifndef GENE
#define GENE
#include <genomestore.pb.h>
#include <indexer.hpp>
#include <set>

/**
 * Goals intervals, genes and store genes by genomes
 * 1, store refgenes into individual genes, serialize/unserialize ids
 * 2, initialize interval from genes, serialize genes and intervals
 * 3, initialize a multiset from gene ids
 **/

struct inv {
  inv() { _inv.set_len(-1); };
  inv(uint32_t start, uint32_t len, bool strand) {
    _inv.set_begin(start);
    _inv.set_len(len);
    _inv.set_strand(strand);
  }
  inv(const std::string &ref, const std::string &chr, uint32_t start,
      uint32_t len, bool strand) {
    inv(start, len, strand);
    _inv.set_ref(ref);
    _inv.set_chr(chr);
  };
  const uint32_t &start() { return _inv.get_start(); };
  const uint32_t &len() { return _inv.get_len(); };
  const uint32_t &end() { return _inv.get_start() + _inv.get_len(); };
  const bool &strand() { return _inv.get_strand(); };
  const std::string &ref() { return _inv.get_ref(); };
  const std::string &chr() { return _inv.get_chr(); };
  const std::string &seqs() { return _inv.get_seqs(); };
  bool operator==(const inv &another) {
    return ((ref() == another.ref()) && (chr() == another.chr()) &&
            (start() == another.start()) && (len() == another.len()) &&
            (strand() == another.strand()));
  };
  bool operator!=(const inv &another) { return !((*this) == (another)); };
  std::string info() {
    char tmp[256];
    auto fm = sprintf(tmp, "[%s, %s]: (%d, %d) @ %s @ {strand: %d}",
                      ref().c_str(), chr().c_str(), start(), end(),
                      seqs().substr(50).c_str(), (int)strand());
    return std::string(tmp);
  };
  bool operator<(const inv &another) {
    if (chr() != another.chr())
      return chr() < another.chr();
    if (another.start() > end())
      return true;
    if (another.end() < start())
      return false;
    return (start() < another.start());
  }
  Interval _inv;
};

// bool combine_intervals(interval &left, interval &right) {
//   if (left.ref().find(right.ref()) == std::string::npos)
//     return false;
//   if (left.chr().find(right.chr()) == std::string::npos)
//     return false;
//   if (left.start > right.start)
//     left.start = right.start;
//   if (left.end < right.end)
//     left.end = right.end;
//   return true;
// }

class gene {
private:
  seqdb *sdb = nullptr;

public:
  std::string sym, id = "";
  unsigned cds_start, cds_end, tx_start, tx_end, idx, chr_;
  bool strand;
  std::multiset<interval, intervalCmp> exons;
  interval inv() { return interval{sdb, chr_, tx_start, tx_end, idx, strand}; };
  interval cds() {
    return interval{sdb, chr_, cds_start, cds_end, idx, strand};
  };
  std::string chr() const { return sdb->chrs[chr_]; };
  std::string ref() const { return sdb->name; };
  unsigned long sz() const { return sdb->sizes[chr()]; };
  std::string info() { return inv().info(); };
  std::vector<interval> get_exons() {
    std::vector<interval> out;
    for (auto &it : exons) {
      out.push_back(it);
    }
    return out;
  };
  std::vector<interval> get_introns() {
    std::vector<interval> out;
    if (exons.size() < 2)
      return out; // no intron in this situation
    auto it = exons.begin();
    auto itt = ++it;
    --it;
    while (itt != exons.end()) {
      out.push_back(
          interval{sdb, chr_, it->end + 1, itt->start - 1, idx, strand});
      ++it;
      ++itt;
    }
    return out;
  };
  interval utr(const bool s) {
    unsigned start, end;
    if (s) {
      start = tx_start;
      end = cds_start - 1;
    } else {
      start = cds_end + 1;
      end = tx_end;
    }
    return interval{sdb, chr_, start, end, idx, strand};
  };
  interval utr5() {
    if (noncoding()) {
      // std::cerr << "Gene at " << idx << " is noncoding, no 5'UTR" <<
      // std::endl;
      return interval();
    }
    return utr(strand);
  };
  interval utr3() {
    if (noncoding()) {
      // std::cerr << "Gene at " << idx << " is noncoding, no 3'UTR" <<
      // std::endl;
      return interval();
    }
    return utr(!strand);
  };
  interval get_p(const int &l, const int &r, const bool &s) {
    unsigned start, end;
    if (s) {
      start = tx_start - l;
      end = tx_start + r;
    } else {
      start = tx_end + l;
      end = tx_end - r;
    }
    return interval{sdb, chr_, start, end, idx, strand};
  };
  interval get_promoter(const int &l, const int &r) {
    if (noncoding()) {
      std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                << std::endl;
      return interval();
    }
    return get_p(l, r, strand);
  };
  interval get_tail(const int &l, const int &r) {
    if (noncoding()) {
      std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                << std::endl;
      return interval();
    }
    return get_p(l, r, !strand);
  };
  virtual bool noncoding() {
    if (cds_start == cds_end) {
      return true;
    }
    return false;
  };
  bool operator==(const gene &another) { return (idx == another.idx); };
  bool operator!=(const gene &another) { return (idx != another.idx); };
};
#endif