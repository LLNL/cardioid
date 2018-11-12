#ifndef STRING_UTILS_HH
#define STRING_UTILS_HH

#include <string>
#include <vector>

/** Concatenates the strings in vv into a single string with a space
 * separating each element.  No space is added to the beginning or
 * end. */
std::string concat(const std::vector<std::string> vv);

#endif
