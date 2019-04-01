#ifndef GENE
#define GENE
#include "interval.hpp"
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
  std::unique_ptr<genomestore::Gene> _gene =
      std::make_unique<genomestore::Gene>();
  std::multiset<inv> introns, exons;

public:
  gene(){};
  gene(std::string &&id) { _gene->set_id(id); };
  gene(std::string &&id, std::string &&sym) : gene(std::move(id)) {
    _gene->set_sym(sym);
  };
  gene(std::string &&id, std::string &&sym, std::string &&ref,
       std::string &&chr)
      : gene(std::move(id), std::move(sym)) {
    _gene->set_ref(ref);
    _gene->set_chr(chr);
  };
  void set_tx(inv &&ii) { _gene->mutable_tx()->CopyFrom(ii.data()); };
  void set_cds(inv &&ii) { _gene->mutable_cds()->CopyFrom(ii.data()); };
  const inv tx() const { return inv{_gene->tx()}; };
  const inv cds() const { return inv{_gene->cds()}; };
  const bool strand() const { return cds().strand(); };
  std::string chr() { return _gene->chr(); };
  std::string ref() { return _gene->ref(); };
  std::string info() { return cds().info(); };
  std::vector<inv> get_exons();
  std::vector<inv> get_introns();
  void add_exon(inv &&i) {
    auto exp = _gene->mutable_exons();
    exp->Add()->CopyFrom(i.data());
  };
  void init_exons() {
    for (auto &i : _gene->exons()) {
      exons.insert(inv{i});
    }
  };
  void add_intron(inv &&i) {
    auto exp = _gene->mutable_introns();
    exp->Add()->CopyFrom(i.data());
  };
  void init_introns() {
    for (auto &i : _gene->introns()) {
      introns.insert(inv{i});
    }
  };
  inv utr(const bool);
  inv utr5();
  inv utr3();
  inv get_p(const int &, const int &, const bool &);
  inv get_promoter(const int &, const int &);
  inv get_tail(const int &, const int &);
  bool noncoding();
  const bool operator==(const gene &another) const {
    return (cds() == another.cds());
  };
  const bool operator!=(const gene &another) const {
    return (cds() != another.cds());
  };
};
#endif