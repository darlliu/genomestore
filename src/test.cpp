#define CATCH_CONFIG_MAIN
#include "base.hpp"
#include "interval.hpp"
#include <catch2/catch.hpp>

TEST_CASE("TESTING BASE DB IO FUNCTIONS") {
  auto bdb = basedb("testdb", "test.db");
  SECTION("TESTING DB INITIALIZATION") {
    REQUIRE(bdb.getldb().cursor != nullptr);
  }
  SECTION("TESTING DB PUT STRING KEY AND VAL") {
    bdb.setdb("testkey", "testval");
    auto msg = bdb.getdb("testkey");
    REQUIRE(msg == "testval");
  }
}

TEST_CASE("TESTING INTERVAL") {
  SECTION("TESTING INTERVAL INITIALIZATION") {
    auto k = inv();
    REQUIRE(k.null());
    REQUIRE(k.empty() == false);
    auto i = inv (1000000, 50, true);
	auto j = inv ( 1000000, 50, true) ;
    REQUIRE( (i.len() == 50));
    REQUIRE( (i == j) == true);
    //inv j = inv{"mm10", "chr1", 1000050, 1000100, true};
    //REQUIRE(i != j);
  }
  //SECTION("TESTING INTERVAL OPEARTIONS") {
  //  inv i = inv{"mm10", "chr1", 1000000, 1000050, true};
  //  inv j = inv{"mm10", "chr1", 1000050, 1000100, true};
  //  inv k = inv{"mm10", "chr1", 1000025, 1000050, true};
  //  REQUIRE(i + j == inv{"mm10", "chr1", 1000000, 1000100, true});
  //  REQUIRE(i - k == inv{"mm10", "chr1", 1000000, 1000025, true});
  //  REQUIRE(i / k == inv{"mm10", "chr1", 1000025, 1000050, true});
  //}
}
