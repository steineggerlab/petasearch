//
// Created by matchy on 2/1/22.
//

#include <cstring>
#include <cstdlib>
#include "SRAUtil.h"
#include "Util.h"

void strrev(char *strRev, const char *str, int len) {
    int start = 0;
    int end = len - 1;
    while (LIKELY(start <= end)) {
        strRev[start] = str[end];
        strRev[end] = str[start];
        ++start;
        --end;
    }
    strRev[len] = '\0';
}

char *substr(char *origStr, unsigned int start, unsigned int end) {
    char *subStr = static_cast<char *>(calloc(end - start + 1, sizeof(char)));
    strncpy(subStr, origStr + start, end - start);
    return subStr;
}

void stripInvalidChars(const char *src, char *dest) {
    size_t j, n = strlen(src);
    for (size_t i = j = 0; i < n; i++) {
        // FIXME: this branch is likely not useful now
        if (src[i] == '\n' || src[i] == '@') {
            continue;
        } else if (src[i] == '*') {
            dest[j++] = 'X';
        } else {
            dest[j++] = src[i];
        }
    }
    dest[j] = '\0';
}

void stripInvalidChars(char *src) {
    size_t j, n = strlen(src);
    for (size_t i = j = 0; i < n; i++) {
        if (src[i] == '\n' || src[i] == '@') {
            continue;
        } else if (src[i] == '*') {
            src[j++] = 'X';
        } else {
            src[j++] = src[i];
        }
    }
    src[j] = '\0';
}
