#define CATCH_CONFIG_MAIN
#include "base.hpp"
#include "gene.hpp"
#include "genomestore.pb.h"
#include "interval.hpp"
#include <catch2/catch.hpp>

TEST_CASE("TESTING BASE DB IO FUNCTIONS") {
  auto bdb = basedb("testdb", "test.db");
  SECTION("TESTING DB INITIALIZATION") {
    const void *ptr = &bdb.getldb();
    REQUIRE(ptr != nullptr);
  }
  SECTION("TESTING DB PUT STRING KEY AND VAL") {
    bdb.setdb("testkey", "testval");
    auto msg = bdb.getdb("testkey");
    REQUIRE(msg == "testval");
  }
}

TEST_CASE("TESTING PROTOS") {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  SECTION("TESTING INTERVAL PROTOS") {
    std::string chr = "chr1", ref = "mm9";
    auto k = genomestore::Interval{};
    REQUIRE(k.len() == 0);
    k.set_len(50);
    REQUIRE(k.len() == 50);
    k.set_start(10000050);
    REQUIRE(k.start() == 10000050);
    k.set_chr(chr);
    REQUIRE((k.chr() == chr));
    k.set_ref(ref);
    REQUIRE((k.chr() == chr));
    REQUIRE((k.ref() == ref));
  }
}

TEST_CASE("TESTING INTERVAL") {
  SECTION("TESTING INTERVAL INITIALIZATION") {
    auto k = inv();
    REQUIRE(k.null());
    REQUIRE(k.empty() == false);
    auto i = inv(1000000, 50, true);
    auto j = inv{1000000, 50, true};
    REQUIRE((i.len() == 50));
    REQUIRE((i == j) == true);
    inv l = inv{"mm10", "chr1", 1000050, 1000100, true};
    REQUIRE((i != l) == true);
  }
  SECTION("TESTING INTERVAL OPEARTIONS") {
    inv i = inv{"mm10", "chr1", 1000000, 1000050, true};
    inv j = inv{"mm10", "chr1", 1000050, 1000100, true};
    inv k = inv{"mm10", "chr1", 1000025, 1000050, true};
    REQUIRE((i + j == inv{"mm10", "chr1", 1000000, 1000100, true}) == true);
    REQUIRE((i - k == inv{"mm10", "chr1", 1000000, 1000025, true}) == true);
    REQUIRE((i / k == inv{}) == true);
    REQUIRE((k / j == j) == true);
  }
  SECTION("TESTING INTERVAL SERIALIZE") {
    auto bdb = basedb("testdb", "test.db");
    inv i = inv{"mm10", "chr1", 1000000, 1000050, true};
    serialize_to_db(bdb, "testinterval", i.data());
    inv j{};
    deserialize_from_db(bdb, "testinterval", j.data());
    REQUIRE((i == j) == true);
  }
  SECTION("TESTING GENE") {
    inv cds = inv{"mm10", "chr1", 1000000, 1000050, true};
    inv tx = inv{"mm10", "chr1", 1000010, 1000040, true};
    inv ex1 = inv{"mm10", "chr1", 1000010, 1000020, true};
    inv ex2 = inv{"mm10", "chr1", 1000030, 1000040, true};
    gene g = gene{"rf1", "testGene", "mm10", "chr1"};
    g.set_cds(inv{"mm10", "chr1", 1000000, 1000050, true});
    g.set_tx(inv{"mm10", "chr1", 1000010, 1000040, true});
    g.add_exon(inv{"mm10", "chr1", 1000010, 1000020, true});
    g.add_exon(inv{"mm10", "chr1", 1000030, 1000040, true});
    g.init_exons();
    REQUIRE((g.cds() == cds) == true);
    REQUIRE((g.tx() == tx) == true);
    REQUIRE((g.get_exons()[0] == ex1) == true);
    REQUIRE((g.get_exons()[1] == ex2) == true);
    REQUIRE(g.get_introns()[0].info() != ex1.info());
  }
}
