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
  auto _gene = std::make_unique<genomestore::Gene>();
  std::multiset<inv,
                [](const inv &left, const inv &right) { return left < right; }>
      introns, exons;

public:
  gene(){};
  inv tx() { return inv{_gene->tx()}; };
  inv cds() { return inv{_gene->cds()}; };
  std::string chr() const { return _gene->chr(); };
  std::string ref() const { return _gene->ref(); };
  std::string info() { return cds().info(); };
  std::vector<inv> get_exons() {
    std::vector<inv> out;
    for (auto &it : exons) {
      out.push_back(it);
    }
    return out;
  };
  std::vector<inv> get_introns() {
    std::vector<inv> out;
    auto exons = get_exons();
    if (exons.size() < 2)
      return out; // no intron in this situation
    auto it = exons.begin();
    auto itt = ++it;
    --it;
    while (itt != exons.end()) {
      out.push_back(inv{ref(), chr(), it->end() + 1, itt->start() - 1, strand});
      ++it;
      ++itt;
    }
    return out;
  };
  inv utr(const bool s) {
    uint32_t start, end;
    if (s) {
      start = tx().start();
      end = cds().start() - 1;
    } else {
      start = cds().end() + 1;
      end = tx().end();
    }
    return inv{ref(), chr(), start, end, strand};
  };
  inv utr5() {
    if (noncoding()) {
      // std::cerr << "Gene at " << idx << " is noncoding, no 5'UTR" <<
      // std::endl;
      return inv{};
    }
    return utr(strand);
  };
  inv utr3() {
    if (noncoding()) {
      // std::cerr << "Gene at " << idx << " is noncoding, no 3'UTR" <<
      // std::endl;
      return inv{};
    }
    return utr(!strand);
  };
  inv get_p(const int &l, const int &r, const bool &s) {
    uint32_t start, end;
    if (s) {
      start = tx().start() - l;
      end = tx().start() + r;
    } else {
      start = tx().end() + l;
      end = tx().end() - r;
    }
    return inv{ref(), chr(), start, end, strand};
  };
  interval get_promoter(const int &l, const int &r) {
    if (noncoding()) {
      std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                << std::endl;
      return inv{};
    }
    return get_p(l, r, strand);
  };
  interval get_tail(const int &l, const int &r) {
    if (noncoding()) {
      std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                << std::endl;
      return inv{};
    }
    return get_p(l, r, !strand);
  };
  virtual bool noncoding() {
    if (cds_start == cds_end) {
      return true;
    }
    return false;
  };
  bool operator==(const gene &another) { return (cds() == another.cds()); };
  bool operator!=(const gene &another) { return (cds() != another.cds()); };
};
#endif