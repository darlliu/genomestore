#ifndef INTERVAL
#define INTERVAL
#include "genomestore.pb.h"
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
    _inv.set_start(start);
    _inv.set_len(len);
    _inv.set_strand(strand);
  }
  inv(const std::string &ref, const std::string &chr, uint32_t start,
      uint32_t end, bool strand) {
    assert(end >= start);
    inv(start, end - start, strand);
    _inv.set_ref(ref);
    _inv.set_chr(chr);
  };
  inv(const std::string &ref, const std::string &chr, uint32_t start,
      uint32_t end, bool strand, const std::string &seqs) {
    inv(ref, chr, start, end, strand);
    _inv.set_seqs(seqs);
  };
  const uint32_t len() const { return _inv.len(); };
  const bool null() const { return len() == -1; };
  const bool empty() const { return len() == 0; };
  const uint32_t start() const { return _inv.start(); };
  const uint32_t end() const { return _inv.start() + _inv.len(); };
  const bool strand() const { return _inv.strand(); };
  const std::string ref() const { return _inv.ref(); };
  const std::string chr() const { return _inv.chr(); };
  const std::string seqs() const { return _inv.seqs(); };
  const std::string info() const;
  const bool operator==(const inv &another);
  const bool operator!=(const inv &another) { return !(*this == another); };
  const bool operator<(const inv &another);
  inv operator+(const inv &another);
  inv operator-(const inv &another);
  inv operator/(inv &another) { return (another - (*this)); };

private:
  genomestore::Interval _inv;
};

#endif