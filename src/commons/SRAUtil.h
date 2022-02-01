//
// Created by matchy on 2/1/22.
//

#ifndef SRASEARCH_SRAUTIL_H

/**
 * @brief Reverse a string and put the results in strRev
 * @param strRev the destination to put the reversed string
 * @param str the string to be revesed
 * @param len the length of the string
 */
void strrev(char *strRev, const char *str, int len);

/**
 * @brief Make a slice origStr[start:end], start inclusive, end exclusive
 * */
char *substr(char *origStr, unsigned int start, unsigned int end);

/**
 * @brief Strip invalid characters from a string (@, *, newline, tab, etc)
 */
void stripInvalidChars(const char *src, char *dest);

#define SRASEARCH_SRAUTIL_H

#endif //SRASEARCH_SRAUTIL_H
