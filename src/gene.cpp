#ifndef GENE_CPP
#define GENE_CPP
#include "gene.hpp"

std::vector<inv> gene::get_exons() {
  std::vector<inv> out;
  for (auto &it : exons) {
    out.push_back(inv{ref(), chr(), it.start(), it.end(), strand()});
  }
  return out;
}

std::vector<inv> gene::get_introns() {
  std::vector<inv> out;
  auto exons = get_exons();
  if (exons.size() < 2)
    return out; // no intron in this situation
  auto it = exons.begin();
  auto itt = ++it;
  --it;
  while (itt != exons.end()) {
    out.push_back(inv{ref(), chr(), it->end() + 1, itt->start() - 1, strand()});
    ++it;
    ++itt;
  }
  return out;
}

inv gene::utr(const bool s) {
  uint32_t start, end;
  if (s) {
    start = tx().start();
    end = cds().start() - 1;
  } else {
    start = cds().end() + 1;
    end = tx().end();
  }
  return inv{ref(), chr(), start, end, strand()};
}

inv gene::utr5() {
  if (noncoding()) {
    // std::cerr << "Gene at " << idx << " is noncoding, no 5'UTR" <<
    // std::endl;
    return inv{};
  }
  return utr(strand());
}

inv gene::utr3() {
  if (noncoding()) {
    // std::cerr << "Gene at " << idx << " is noncoding, no 3'UTR" <<
    // std::endl;
    return inv{};
  }
  return utr(!strand());
}

inv gene::get_p(const int &l, const int &r, const bool &s) {
  uint32_t start, end;
  if (s) {
    start = tx().start() - l;
    end = tx().start() + r;
  } else {
    start = tx().end() + l;
    end = tx().end() - r;
  }
  return inv{ref(), chr(), start, end, strand()};
}

inv gene::get_promoter(const int &l, const int &r) {
  if (noncoding()) {
    return inv{};
  }
  return get_p(l, r, strand());
}

inv gene::get_tail(const int &l, const int &r) {
  if (noncoding()) {
    return inv{};
  }
  return get_p(l, r, !strand());
}

bool gene::noncoding() {
  if (cds().start() == cds().end()) {
    return true;
  }
  return false;
}
#endif