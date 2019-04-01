#ifndef INTERVAL
#define INTERVAL
#include "base.hpp"
#include "genomestore.pb.h"
#include <set>

/**
 * Goals intervals, genes and store genes by genomes
 * 1, store refgenes into individual genes, serialize/unserialize ids
 * 2, initialize interval from genes, serialize genes and intervals
 * 3, initialize a multiset from gene ids
 **/

struct inv {
  inv() { _inv->set_len(-1); };
  inv(genomestore::Interval ii) { *_inv = ii; };
  inv(uint32_t start, uint32_t len, bool strand) {
    _inv->set_start(start);
    _inv->set_len(len);
    _inv->set_strand(strand);
  }
  inv(std::string &&ref, std::string &&chr, uint32_t start, uint32_t end,
      bool strand)
      : inv(start, end - start, strand) {
    assert(end >= start);
    auto _ref = ref, _chr = chr;
    _inv->set_ref(std::move(_ref));
    _inv->set_chr(std::move(_chr));
  };
  inv(std::string &&ref, std::string &&chr, uint32_t start, uint32_t end,
      bool strand, std::string &&seqs)
      : inv(std::move(ref), std::move(chr), start, end, strand) {
    _inv->set_seqs(std::move(seqs));
  };
  genomestore::Interval &data() { return *_inv; };
  const uint32_t len() const { return _inv->len(); };
  const bool null() const { return len() == -1; };
  const bool empty() const { return len() == 0; };
  const uint32_t start() const { return _inv->start(); };
  const uint32_t end() const { return _inv->start() + _inv->len(); };
  const bool strand() const { return _inv->strand(); };
  const std::string ref() const { return _inv->ref(); };
  const std::string chr() const { return _inv->chr(); };
  const std::string seqs() const { return _inv->seqs(); };
  const std::string info() const;
  const bool operator==(const inv &another) const;
  const bool operator!=(const inv &another) const {
    return !(*this == another);
  };
  const bool operator<(const inv &another) const;
  inv operator+(const inv &another);
  inv operator-(const inv &another);
  inv operator/(inv &another) { return (another - (*this)); };

private:
  std::unique_ptr<genomestore::Interval> _inv =
      std::make_unique<genomestore::Interval>();
};
#endif