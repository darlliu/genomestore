#ifndef INTERVAL_CPP
#define INTERVAL_CPP
#include "interval.hpp"

const bool inv::operator==(const inv &another) {
  return ((ref() == another.ref()) && (chr() == another.chr()) &&
          (start() == another.start()) && (len() == another.len()) &&
          (strand() == another.strand()));
}

const std::string inv::info() const {
  char tmp[256];
  auto fm = sprintf_s(tmp, "[%s, %s]: (%d, %d) @ %s @ {strand: %d}",
                      ref().c_str(), chr().c_str(), start(), end(),
                      seqs().substr(50).c_str(), (int)strand());
  return std::string(tmp);
}

const bool inv::operator<(const inv &another) {
  if (chr() != another.chr())
    return chr() < another.chr();
  if (another.start() > end())
    return true;
  if (another.end() < start())
    return false;
  return (start() < another.start());
}

inv inv::operator+(const inv &right) {
	auto _ref = ref(), _chr = chr();
  // combine intervals, while taking the strandity of the current interval
  if (ref().find(right.ref()) == std::string::npos ||
      chr().find(right.chr()) == std::string::npos) {
	  return inv{ std::move(_ref), std::move(_chr), start(), end(), strand()};
  }
  uint32_t _start, _end;
  if (start() > right.start())
    _start = right.start();
  else
    _start = start();
  if (end() < right.end())
    _end = right.end();
  else
    _end = end();

  return inv{std::move(_ref), std::move(_chr), _start, _end, strand()};
}

inv inv::operator-(const inv &right) {
	auto _ref = ref(), _chr = chr();
	// combine intervals, while taking the strandity of the current interval
	if (ref().find(right.ref()) == std::string::npos ||
		chr().find(right.chr()) == std::string::npos) {
		return inv{ std::move(_ref), std::move(_chr), start(), end(), strand() };
	}
  uint32_t _start, _end;
  if (start() > right.start())
    _start = start();
  else
    _start = right.start();
  if (end() < right.end())
    _end = end();
  else
    _end = right.end();
  return inv{std::move(_ref), std::move(_chr), _start, _end, strand()};
}
#endif