/*
 * A base indexer for sequentially indexed contents
 * contains a map from string (by chromosome for example) to vector of indices
 * and a sequentially sorted vector with a fixed step size from indices to keys
 * and a on disk hash table from keys to actual content (chunks of sequences)
 * Functionalities are provided to query sequences from any slice of indices on
 * major key
 */

#ifndef INDEXER
#define INDEXER
#include "base.hpp"
#include <sstream>

class seqdb {

public:
  seqdb(const std::string &name, const uint32_t &sz,
        const std::string &dbp = "./test")
      : name(name), chunksz(sz), dbpath(dbp), chr(""), assemble(false),
        scaffold(false){};
  seqdb() : seqdb("default", 10000){};
  ~seqdb();
  virtual void import(const std::string &);
  virtual void scan(const std::string &);
  virtual void import_scan(const std::string &);
  virtual void import_feed();
  // import from a fasta file and build a db
  virtual void import_chr_fsq();
  // import a chromosome file
  virtual void import_chr(const std::string &);
  std::string get3(const std::string &chr, const uint32_t &l,
                   const uint32_t &r) {
    return this->get(chr, l, r); // this is an adaptor for python extension
  };
  virtual std::string get(const std::string &chr, const uint32_t &l,
                          const uint32_t &r);
  virtual std::string get(const uint32_t &, const uint32_t &);
  virtual std::string get(const std::string &key) { return ""; };
  virtual std::string get(const interval &itv) {
    return get(itv.chr, itv.l, itv.r);
  };
  void close_db(DB &dbs) { dbs.clear(); }; // nothing needs to be done!
  virtual void load_db(const std::string &);
  virtual void load_db_();
  void load_sizes(std::ifstream &);
  virtual void export_db(const std::string &);
  void set_chr(const std::string &c) { chr = c; };
  void set(const std::string &, const std::string &);
  void del(const std::string &);
  virtual uint32_t get_index(const uint32_t &idx) {
    return idx / chunksz * chunksz;
  };
  void serialize(const std::string &fname) { export_db(fname); };
  uint32_t chunk_sz() { return chunksz; };

  std::string name,
      dbpath;       // name of the db and the directory path for kyotocabinet
  uint32_t chunksz; // size of each chunk of sequence
  std::string chr;
  std::vector<std::string> chrs;
  INDEXMAP indices;
  DB dbs;
  SIZES sizes, shifts, postfixes;
  NAMES fapaths, dbpaths;
  bool assemble, scaffold;
  std::vector<std::string> tmp_dbpaths;
  DBFSQ dbs_fsq, dbs_fsz;
};

#endif
