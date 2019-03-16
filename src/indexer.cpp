#include "indexer.hpp"
using namespace std::experimental::filesystem;

seqdb::~seqdb() {
  close_db(dbs);
  for (auto &it : dbs_fsq) {
    it.second->close();
  }
}

std::string seqdb::get(const std::string &chr, const uint32_t &l,
                       const uint32_t &r) {
  if (sizes.count(chr) == 0 || sizes[chr] <= r || r <= l)
    return "";
  set_chr(chr);
  return get(l, r);
}

void seqdb::scan(const std::string &dirname) {
  auto dp = path(dirname);
  chrs.clear();
  std::vector<path> fps;
  try {
    if (exists(dp) && is_directory(dp)) {
      std::cerr << " Now reading " << dirname << std::endl;
    } else
      throw("Path does not exist: " + dirname);
  } catch (const filesystem_error &ex) {
    std::cerr << "Error opening path " << dirname << std::endl;
    throw(ex);
  }
  for (auto it = directory_iterator(dp); it != directory_iterator(); ++it) {
    if (!is_directory(it->path())) {
      auto fp = it->path();
      if (fp.extension().string().find("fa") == std::string::npos)
        continue;
      std::cout << "Found fasta file : " << fp.filename() << std::endl;
      fps.push_back(fp);
    } else {
      std::cout << "skipping directory : " << it->path().filename()
                << std::endl;
    }
  }
  std::cout << "Found files count: " << fps.size() << std::endl;
  for (unsigned i = 0; i < fps.size(); ++i) {
    auto fp = fps[i];
    auto chr = fp.filename().stem().string();
    chrs.push_back(chr);
    fapaths[chr] = fp.string();
    auto dbp = dbpath + "/" + chr + ".fsq";
    dbpaths[chr] = dbp;
  }
}

void seqdb::import(const std::string &dirname) {
  scan(dirname);
  for (auto &chr : chrs) {
    set_chr(chr);
    import_chr_fsq();
  }
  load_db_();
}

void seqdb::import_scan(const std::string &dirname) {
  scan(dirname);
  load_db_();
}

void seqdb::import_chr(const std::string &fname) {
  auto fp = path(fname);
  if (exists(fp) && fp.extension().string().find("fa") != std::string::npos) {
    std::cerr << " Now reading " << fname << std::endl;
  } else
    throw("File does not exist or wrong type: " + fname);
  fapaths[chr] = fp.string();
  auto chr = fp.filename().stem().string();
  chrs.push_back(chr);
  set_chr(chr);
  auto dbp = dbpath + "/" + chr + ".fsq";
  import_chr_fsq();
  return;
}
void seqdb::import_chr_fsq() {
  auto chrfp = path(fapaths[chr]);
  auto dbp = path(dbpaths[chr]);
  std::ofstream dbf1(dbp, std::ofstream::binary);
  dbp.replace_extension("szi");
  std::ofstream dbf2(dbp, std::ofstream::binary);
  std::ifstream ifs(chrfp);
  if (!ifs.is_open() || !dbf1.is_open() || !dbf2.is_open())
    throw("Error Opening file import sequence!");
  uint32_t idx = 0;
  std::string tmp = "", line;
  std::cout << "Opened file " << chrfp << " On chr: " << chr << std::endl;
  while (std::getline(ifs, line)) {
    if (line.size() == 0 || line[0] == '>')
      continue;
    tmp += line;
  }
  ifs.close();
  auto buf = new char[(tmp.size() + 1) / 2];
  encode_seq(tmp, buf);
  auto pos = dbf1.tellp();
  assert(pos < UINT32_MAX);
  dbf1.write(buf, (tmp.size() + 1) / 2);
  if ((!assemble) && !(dbf1.good()))
    throw("Error writing to a 4bit sequence " + chr);
  dbf1.close();
  assert(tmp.size() < UINT32_MAX);
  idx += static_cast<uint32_t>(tmp.size());
  sizes[chr] = idx;
  shifts[chr] = static_cast<uint32_t>(pos);
  dbf2 << (chr + " " + std::to_string(pos) + " " + std::to_string(idx) + " ");
  std::cout << "Loaded length " << idx << std::endl;
  dbf2.close();
  delete[] buf;
  return;
};

void seqdb::import_feed() {
  std::cerr << "..... Now importing from stdin ....." << std::endl;
  std::string line, tmp("");
  std::ofstream dbf1, dbf2;
  path dbp;
  if (scaffold) {
    dbp = path(dbpath + "/" + name + ".fsq");
    dbf1.open(dbp, std::ofstream::binary);
    dbpaths[name] = dbp.string();
    dbp.replace_extension("szi");
    dbf2.open(dbp, std::ofstream::binary);
  }
  uint32_t idx = 0;
  bool flag = false;
  auto inner = [&]() {
    auto buf = new char[(tmp.size() + 1) / 2];
    encode_seq(tmp, buf);
    auto pos = dbf1.tellp();
    dbf1.write(buf, (tmp.size() + 1) / 2);
    assert(tmp.size() < UINT32_MAX);
    idx += static_cast<uint32_t>(tmp.size());
    dbf2 << (chr + " " + std::to_string(pos) + " " + std::to_string(idx) + " ");
    std::cout << "Loaded length " << idx << std::endl;
    idx = 0;
    dbf2.flush();
    tmp = "";
    delete[] buf;
  };
  while (std::cin) {
    std::getline(std::cin, line);
    if (line.size() == 0)
      continue;
    if (line[0] == '>') {
      if (flag)
        inner();
      chr = line.substr(1);
      std::cerr << "Found breaking point " << chr
                << " . Assuming this is a chromosome-like!" << std::endl;
      if (!scaffold) {
        dbf1.close();
        dbf2.close();
        dbp = path(dbpath + "/" + chr + ".fsq");
        dbpaths[chr] = dbp.string();
        dbf1.open(dbp, std::ofstream::binary);
        dbp.replace_extension("szi");
        dbf2.open(dbp, std::ofstream::binary);
      }
      flag = true;
      continue;
    }
    tmp += line;
  }
  inner();
  dbf1.close();
  dbf2.close();
  load_db_();
  std::cerr << "Finished importing" << std::endl;
  return;
}

std::string seqdb::get(const uint32_t &l, const uint32_t &r) {
  if (r < l)
    throw("Interval incorrect!");
  std::shared_ptr<std::ifstream> db;
  unsigned pos = 0;
  if (scaffold) {
    db = dbs_fsq[name];
    pos = shifts[chr];
  } else
    db = dbs_fsq[chr];
  unsigned ll, rr;
  if (l % 2)
    ll = l - 1;
  else
    ll = l;
  if (r % 2)
    rr = r + 1;
  else
    rr = r;
  unsigned rr2 = rr / 2, ll2 = ll / 2;
  db->seekg(0, db->beg); // seek back
  db->seekg(pos + ll2, db->beg);
  auto buf = new char[rr2 - ll2];
  db->read(buf, rr2 - ll2);
  db->sync();
  std::string decoded_seq = decode_seq(buf, rr2 - ll2);
  delete[] buf;
  return decoded_seq.substr(l - ll, r - l);
};

void seqdb::export_db(const std::string &fp) {
  std::cerr << "Trying to serialize into " << fp << std::endl;
  std::ofstream ofs(fp);
  if (!ofs.is_open())
    throw("open error (export db): " + fp);
  json j;
  j["name"] = json(name);
  j["dbpath"] = json(dbpath);
  j["chunksz"] = json(chunksz);
  j["chrs"] = json(chrs);
  j["indices"] = json(indices);
  j["sizes"] = json(sizes);
  j["shifts"] = json(shifts);
  j["postfixes"] = json(postfixes);
  j["faptahs"] = json(fapaths);
  j["dbpaths"] = json(dbpaths);
  j["scaffold"] = json(scaffold);
  j["assemble"] = json(assemble);
  ofs << j;
  std::cerr << "Writing out SeqDB ... " << std::endl;
};
void seqdb::load_sizes(std::ifstream &ifs) {
  std::string ch;
  unsigned pos;
  std::stringstream ss;
  ss << ifs.rdbuf();
  while (ss.good()) {
    ss >> ch;
    if (!ss.good())
      break;
    ss >> pos;
    shifts[ch] = pos;
    ss >> pos;
    sizes[ch] = pos;
    chrs.push_back(ch);
  }
  ifs.close();
  return;
}

void seqdb::load_db_() {
  if (scaffold) {
    dbs_fsq[name] = std::shared_ptr<std::ifstream>(
        new std::ifstream(dbpath + "/" + name + ".fsq",
                          std::ifstream::in | std::ifstream::binary));
    dbs_fsz[name] = std::shared_ptr<std::ifstream>(
        new std::ifstream(dbpath + "/" + name + ".szi"));
    if (!dbs_fsq[name]->good() || !dbs_fsz[name]->good())
      throw("Error opening a 4bit file at " + name);
    load_sizes(*dbs_fsz[name]);
  } else {
    for (auto it : dbpaths) {
      auto dbp = path(it.second);
      dbs_fsq[it.first] =
          std::shared_ptr<std::ifstream>(new std::ifstream(dbp));
      dbp.replace_extension("szi");
      dbs_fsz[it.first] =
          std::shared_ptr<std::ifstream>(new std::ifstream(dbp));
      if (!dbs_fsq[it.first]->good() || !dbs_fsz[it.first]->good())
        throw("Error opening a 4bit file at " + it.first);
      load_sizes(*dbs_fsz[it.first]);
    }
  }
  return;
};

void seqdb::load_db(const std::string &fp) {
  std::cerr << "Trying to deserialize from " << fp << " ... ";
  std::ifstream ifs(fp);
  if (!ifs.is_open())
    throw("open error (load db): " + fp);
  json j;
  ifs >> j;
  name = j["name"].get<std::string>();
  dbpath = j["dbpath"].get<std::string>();
  chunksz = j["chunksz"].get<uint32_t>();
  chrs = j["chrs"].get<std::vector<std::string>>();
  indices = j["indices"].get<INDEXMAP>();
  sizes = j["sizes"].get<SIZES>();
  shifts = j["shifts"].get<SIZES>();
  postfixes = j["postfixes"].get<SIZES>();
  fapaths = j["faptahs"].get<NAMES>();
  dbpaths = j["dbpaths"].get<NAMES>();
  scaffold = j["scaffold"].get<bool>();
  assemble = j["assemble"].get<bool>();
  std::cerr << "Loading SeqDB ... " << std::endl;
  load_db_();
};