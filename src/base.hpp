#ifndef BASE
#define BASE

#define DEBUG 0

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <lmdb.h>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#if defined(WIN32) || defined(_WIN32) ||                                       \
    defined(__WIN32) && !defined(__CYGWIN__)
#include <experimental/filesystem>
#else
#include <experimental/filesystem>
#endif
using json = nlohmann::json;

typedef std::vector<uint32_t> INDICES;
typedef std::unordered_map<std::string, INDICES> INDEXMAP;
typedef std::unordered_map<std::string, std::string> NAMES;
typedef std::unordered_map<std::string, uint32_t> SIZES;

typedef enum {
  increment = 0,
  custom
} INDEXTYPE; // these are optional identifiers left unused for now
typedef enum { reference = 0, matching, other } INTERVALTYPE;

struct interval {
public:
  std::string ref, chr, seq;
  unsigned l = -1, r = -1;
  INTERVALTYPE tt = reference;
  float score = 0.0;
  bool strand = true;
};

struct ldb {
  MDB_env *env;
  MDB_dbi dbi;
  MDB_txn *txn;
  MDB_stat mst;
  MDB_cursor *cursor, *cur2;
  MDB_cursor_op op;
};

char encode_char(char *in);
void decode_char(const char &in, char *out);
void encode_seq(std::string in, char *out);

std::string decode_seq(char *in, unsigned size);

std::string print_interval(const interval &in);

bool combine_intervals(interval &left, interval &right);

std::string get_reverse_comp(const std::string &in);

class basedb {
  // LDB IO class for sequence db, alignment db, annotation db, and SNP db

public:
  basedb(const std::string &name, const std::string &dbpath);
  ~basedb();
  void dbinit(ldb &dbenv,
              const std::string &dbpath); // should be called by all derived
                                          // classes, will use dbpath, name
  void setdb(ldb &dbenv, const std::string &key, const std::string &val,
             bool append, int &cnt);
  std::string getdb(ldb &dbenv, const std::string &key);
  void dbinit() { dbinit(dbenv, dbpath + "_ldb"); }; // default to one db
  const ldb getldb() const { return dbenv; };
  void setdb(const std::string &key, const std::string &val,
             bool append = false) {
    setdb(dbenv, key, val, append, cnt);
  };
  std::string getdb(const std::string &key) { return getdb(dbenv, key); };

private:
  std::string name, dbpath;
  // names of the db, the directory path for the underlying datastore
  ldb dbenv; // default db environment
  int cnt = 0;
};

typedef MDB_val DBVAL;
typedef std::unordered_map<std::string, std::vector<std::shared_ptr<basedb>>>
    DB;
typedef std::unordered_map<std::string, std::shared_ptr<std::ifstream>> DBFSQ;
#endif
