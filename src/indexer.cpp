#include "indexer.hpp"

void seqdb::scan(const std::string &dirname) {
  using namespace std::experimental::filesystem;
  auto fp = path(dirname);
  chrs.clear();
  std::vector<std::string> fps, fns;
  try {
    if (exists(fp) && is_directory(fp)) {
      std::cerr << " Now reading " << dirname << std::endl;
    } else
      throw("Path does not exist: " + dirname);
  } catch (const filesystem_error &ex) {
    std::cerr << "Error opening path " << dirname << std::endl;
    throw(ex);
  }
  for (auto it = directory_iterator(fp); it != directory_iterator(); ++it) {

    if (!is_directory(it->path())) {
      auto fp = it->path().string();
      if (fp.find(".fa") == std::string::npos)
        continue;
      std::cout << "Found fasta file : " << fp << std::endl;
      fps.push_back(fp);
      fns.push_back(it->path().filename().string());

    } else {
      std::cout << "skipping directory : " << it->path().filename()
                << std::endl;
    }
  }
  std::cout << "Found files count: " << fps.size() << std::endl;
  for (std::string fn : fns) {
    replace_last(fn, ".fasta", "");
    replace_last(fn, ".fa", "");
    chrs.push_back(fn);
  }
  for (unsigned i = 0; i < fps.size(); ++i) {
    auto chr = chrs[i];
    auto chrp = fps[i];
    fapaths[chr] = chrp;
    auto dbp = dbpath + "/" + chr + ".fsq";
    dbpaths[chr] = dbp;
  }
}

void seqdb::import(const std::string &dirname) {
  scan(dirname);
  for (unsigned i = 0; i < chrs.size(); ++i) {
    auto chr = chrs[i];
    set_chr(chr);
    import_chr_fsq();
  }
  load_db_();
}

void seqdb::import_scan(const std::string& dirname){
  scan(dirname);
  load_db_();
}

void seqdb::import_chr(const std::string &fname) {
  using namespace std::experimental::filesystem;
  auto fp = path(fname);
  if (exists(fp) && fname.find(".fa") != std::string::npos) {
    std::cerr << " Now reading " << fname << std::endl;
  } else
    throw("File does not exist: " + fname);
  auto fn = fp.filename().string();
  fapaths[chr] = fp.string();
  replace_last(fn, ".fasta", "");
  replace_last(fn, ".fa", "");
  chrs.push_back(fn);
  set_chr(fn);
  fapaths[chr] = fp.string();
  auto dbp = dbpath + "/" + chr + ".fsq";
  import_chr_fsq();
  return;
}
void seqdb::import_chr_fsq() {
  auto chrfp = fapaths[chr];
  auto dbp = dbpaths[chr];
  std::ofstream dbf1(dbp, std::ofstream::binary);
  replace_last(dbp, ".fsq", ".szi");
  std::ofstream dbf2(dbp, std::ofstream::binary);
  std::ifstream ifs(chrfp);

  if (!ifs.is_open() || !dbf1.is_open() || !dbf2.is_open())
    throw("Error Opening file import sequence!");
  size_t idx = 0;
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
  dbf1.write(buf, (tmp.size() + 1) / 2);
  if ((!assemble) && !(dbf1.good()))
    throw("Error writing to a 4bit sequence " + chr);
  dbf1.close();
  idx += tmp.size();
  sizes[chr] = idx;
  shifts[chr] = pos;
  dbf2 << (chr + " " + std::to_string(pos) + " " + std::to_string(idx) + " ");
  std::cout << "Loaded length " << idx << std::endl;
  dbf2.close();
  delete[] buf;
  return;
};

void seqdb::import_feed() {
  std::cerr << "..... Now importing from stdin ....." << std::endl;
  std::string line, tmp(""), dbp;
  std::ofstream dbf1, dbf2;
  if (scaffold) {
    dbp = dbpath + "/" + name + ".fsq";
    dbf1.open(dbp, std::ofstream::binary);
    dbpaths[name] = dbp;
    replace_last(dbp, ".fsq", ".szi");
    dbf2.open(dbp, std::ofstream::binary);
  }
  size_t idx = 0;
  bool flag = false;
  auto inner = [&]() {
    auto buf = new char[(tmp.size() + 1) / 2];
    encode_seq(tmp, buf);
    auto pos = dbf1.tellp();
    dbf1.write(buf, (tmp.size() + 1) / 2);
    idx += tmp.size();
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
        dbp = dbpath + "/" + chr + ".fsq";
        dbpaths[chr] = dbp;
        dbf1.open(dbp, std::ofstream::binary);
        replace_last(dbp, ".fsq", ".szi");
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

void seqdb::init_db(const std::vector<std::string> &chrs, DB &dbs) {
  // if (dbs.size() > 0)
  //   close_db(dbs);
  // std::cout << "Initializing DBs, length: " << chrs.size() << std::endl;
  // for (auto chr : chrs) {
  //   dbs[chr].push_back(std::shared_ptr<basedb>(new basedb()));
  // }
  return;
};

std::string seqdb::get(const size_t &l, const size_t &r) {
  if (r < l)
    throw("Interval incorrect!");
  std::shared_ptr<std::ifstream> db;
  unsigned pos = 0;
  // std::cerr <<"Trying to get "<<scaffold <<" "<<l <<", "<<r<<", "<<pos <<",
  // "<<sizes[chr]<<std::endl;
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
  // std::cerr <<"Reading file for get "<<ll2<<", "<<rr2<<" from
  // "<<pos<<"total"<<rr2-ll2<<std::endl;
  db->read(buf, rr2 - ll2);
  db->sync();
  // std::cerr <<"Decoding string for get count "<<db->gcount()<<std::endl;
  std::string decoded_seq = decode_seq(buf, rr2 - ll2);
  delete[] buf;
  // std::cerr <<"returning substr from" <<decoded_seq<<std::endl;
  return decoded_seq.substr(l - ll, r - l);
};

void seqdb::export_db_kch(const std::string &kdbname) {
  export_db(kdbname);
}

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
  ofs << j ;
  std::cerr << "Writing out SeqDB ... " << std::endl;

};
void seqdb::load_sizes(std::ifstream &ifs) {
  std::string ch;
  unsigned pos;
  std::stringstream ss;
  ss << ifs.rdbuf();
  while (ss.good()) {
    ss >> ch;
    if(!ss.good()) break;
    ss >> pos;
    shifts[ch] = pos;
    ss >> pos;
    sizes[ch] = pos;
    chrs.push_back(ch);
    //std::cerr << "Loaded sizes:"<<ch <<", "<<shifts[ch]<<", "<<sizes[ch]<<std::endl;
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
      auto dbp = it.second;
      dbs_fsq[it.first] =
          std::shared_ptr<std::ifstream>(new std::ifstream(dbp));
      replace_last(dbp, ".fsq", ".szi");
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
  chunksz = j["chunksz"].get<size_t>();
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

void seqdb::load_db_kch(const std::string &kdbname, const std::string &key) {
  // TBI
};
