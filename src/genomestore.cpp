#include <indexer.hpp>
#include <pybind11/pybind11.h>
namespace py = pybind11;
PYBIND11_MODULE(genomestore, m) {
  py::class_<seqdb>(m, "seqdb")
      .def(py::init<const std::string &, const size_t, const std::string &>())
      .def("load_genomes", &seqdb::import)
      .def("scan_genomes", &seqdb::import_scan)
      .def("load_genome_stdin", &seqdb::import_feed)
      .def("load_genome", &seqdb::import_chr)
      .def("get", &seqdb::get3)
      .def("load", &seqdb::load_db)
      .def("save", &seqdb::serialize);
}