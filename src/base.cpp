#define DEBUG 0

#include "base.hpp"
#define E(expr) CHECKLMDB((rc = (expr)) == MDB_SUCCESS, #expr)
#define RES(err, expr) ((rc = expr) == (err) || (CHECKLMDB(!rc, #expr), 0))
#define CHECKLMDB(test, msg)                                                   \
  ((test) ? (void)0                                                            \
          : ((void)fprintf(stderr, "%s:%d: %s: %s\n", __FILE__, __LINE__, msg, \
                           mdb_strerror(rc)),                                  \
             abort()))

char encode_char(char *in) {
  std::bitset<8> bs;
  // std::cerr <<"Encoding "<<bs << " from "<<in[0]<<in[1]<<std::endl;
  for (int i = 0; i < 2; ++i)
    switch (in[i]) {
    case 'a': {
      bs[i * 4 + 1] = 1;
    }
    case 'A': {
      bs[i * 4 + 2] = 0;
      bs[i * 4 + 3] = 0;
      break;
    }
    case 't': {
      bs[i * 4 + 1] = 1;
    }
    case 'T': {
      bs[i * 4 + 2] = 1;
      bs[i * 4 + 3] = 1;
      break;
    }
    case 'c': {
      bs[i * 4 + 1] = 1;
    }
    case 'C': {
      bs[i * 4 + 2] = 1;
      bs[i * 4 + 3] = 0;
      break;
    }
    case 'g': {
      bs[i * 4 + 1] = 1;
    }
    case 'G': {
      bs[i * 4 + 2] = 0;
      bs[i * 4 + 3] = 1;
      break;
    }
    default:
      bs[i * 4] = 1;
    }
  // std::cerr <<"Got "<<bs << " and "<<bs.to_ulong()<<std::endl;
  return static_cast<unsigned char>(bs.to_ulong());
}

void decode_char(const char &in, char *out) {
  std::bitset<8> bs((unsigned char)in);
  // std::cerr <<"Decoding "<<bs<<"...";
  for (int i = 0; i < 2; ++i) {
    if (bs[i * 4]) {
      out[i] = 'N';
    } else if (bs[i * 4 + 2] && !bs[i * 4 + 3]) {
      if (bs[i * 4 + 1])
        out[i] = 'c';
      else
        out[i] = 'C';
    } else if (!bs[i * 4 + 2] && bs[i * 4 + 3]) {
      if (bs[i * 4 + 1])
        out[i] = 'g';
      else
        out[i] = 'G';
    } else if (!bs[i * 4 + 2] && !bs[i * 4 + 3]) {
      if (bs[i * 4 + 1])
        out[i] = 'a';
      else
        out[i] = 'A';
    } else if (bs[i * 4 + 2] && bs[i * 4 + 3]) {
      if (bs[i * 4 + 1])
        out[i] = 't';
      else
        out[i] = 'T';
    } else
      throw("Not all cases covered in decoding");
    // std::cerr <<out[i];
  }
  // std::cerr <<std::endl;
  return;
}

void encode_seq(std::string in, char *out) {
  // out.reserve(in.size()/2);
  char tmp[2];
  for (int i = 0; i < (in.size() + 1) / 2; ++i) {
    tmp[0] = in[i * 2];
    if (i * 2 + 1 < in.size())
      tmp[1] = in[i * 2 + 1];
    char c = encode_char(tmp);
    out[i] = c;
  }
  return;
}

std::string decode_seq(char *in, unsigned size) {
  std::string out;
  out.reserve(size * 2);
  for (int i = 0; i < size; ++i) {
    char tmp[2];
    decode_char(in[i], tmp);
    out.push_back(tmp[0]);
    out.push_back(tmp[1]);
  }
  return out;
}

std::string print_interval(const interval &in) {
  char out[256];
  auto seq = in.seq;
  if (seq.size() > 50)
    seq = seq.substr(0, 50) + "..";
  auto fm = sprintf(out, "[%s, %s]: (%d, %d) @ %s @ {score: %f , strand: %d}",
                    in.ref.c_str(), in.chr.c_str(), in.l, in.r, seq.c_str(),
                    in.score, (int)in.strand);
  return std::string(out);
}

bool combine_intervals(interval &left, interval &right) {
  if (left.ref.find(right.ref) == std::string::npos)
    return false;
  if (left.chr.find(right.chr) == std::string::npos)
    return false;
  left.seq = ""; // remove the sequence
  if (left.l > right.l)
    left.l = right.l;
  if (left.r < right.r)
    left.r = right.r;
  return true;
}

std::string get_reverse_comp(const std::string &in) {
  std::string out;
  out.reserve(in.size());
  for (auto &c : in) {
    switch (c) {
    case 'A':
      out.push_back('T');
    case 'T':
      out.push_back('A');
    case 'U':
      out.push_back('A');
    case 'G':
      out.push_back('C');
    case 'C':
      out.push_back('G');
    case 'a':
      out.push_back('t');
    case 't':
      out.push_back('a');
    case 'u':
      out.push_back('a');
    case 'g':
      out.push_back('c');
    case 'c':
      out.push_back('g');
    default:
      out.push_back(c);
    }
  }
  return out;
};

basedb::basedb(const std::string &name, const std::string &dbpath)
    : name(name), dbpath(dbpath) {
  dbinit();
}

basedb::~basedb() {
  // do a few cleaning ops, don't care if fails
  mdb_txn_commit(dbenv.txn);
  mdb_txn_abort(dbenv.txn);
  mdb_dbi_close(dbenv.env, dbenv.dbi);
}

void basedb::dbinit(ldb &dbenv, const std::string &dbpath) {
  int rc;
  E(mdb_env_create(&dbenv.env));
  E(mdb_env_set_maxreaders(dbenv.env, 1));
  E(mdb_env_set_mapsize(dbenv.env, 2 << 18));
  E(mdb_env_open(dbenv.env, dbpath.c_str(), MDB_NOSUBDIR /*|MDB_NOSYNC*/,
                 0664));
  E(mdb_txn_begin(dbenv.env, NULL, 0, &dbenv.txn));
  E(mdb_dbi_open(dbenv.txn, NULL, 0, &dbenv.dbi));
  // opens the db connection as in mtest.c
  return;
}

void basedb::setdb(ldb &dbenv, const std::string &key, const std::string &val,
                   bool append, int &cnt) {
  int rc;
  MDB_val mkey, mdata;
  mkey.mv_size = key.size() * sizeof(char);
  mkey.mv_data = (void *)key.data();

  mdata.mv_size = val.size() * sizeof(char);
  mdata.mv_data = (void *)val.data();

  if (append)
    E(mdb_put(dbenv.txn, dbenv.dbi, &mkey, &mdata, 0 | MDB_APPENDDUP));
  else
    E(mdb_put(dbenv.txn, dbenv.dbi, &mkey, &mdata, 0));
  if (++cnt > 1000) {
    E(mdb_txn_commit(dbenv.txn));
    cnt = 0; // forces a sync every 1000 writes
  };
  return;
}

std::string basedb::getdb(ldb &dbenv, const std::string &key) {
  int rc;
  MDB_val mkey, mdata;
  mkey.mv_size = key.size() * sizeof(char);
  mkey.mv_data = (void *)key.data();

  E(mdb_get(dbenv.txn, dbenv.dbi, &mkey, &mdata));
  int sz = mdata.mv_size / sizeof(char);
  std::string out((char *)mdata.mv_data, sz);
  return out;
}
