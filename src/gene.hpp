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

class gene {
private:
  Gene _gene;
  std::multiset<interval, [](const Interval &left,
                             const Interval &right) { return left < right; }>
      introns, exons;

public:
  gene(){};
  interval tx() { return _gene.tx(); };
  interval cds() { return _gene.cds(); };
  std::string chr() const { return sdb->chrs[chr_]; };
  std::string ref() const { return sdb->name; };
  unsigned long sz() const { return sdb->sizes[chr()]; };
  std::string info() { return tx().info(); };
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