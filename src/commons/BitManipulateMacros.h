//
// Created by matchy233 on 7/3/20.
//

#ifndef SRASEARCH_BITMANIPULATEMACROS_H
#define SRASEARCH_BITMANIPULATEMACROS_H

#define SET_END_FLAG(num)           (0x8000U | (num))
#define IS_LAST_15_BITS(num)        (0x8000U & (num))
#define GET_15_BITS(num)            (0x7fffU & (num))
#define DECODE_15_BITS(diff, num)   ((diff) | (GET_15_BITS(num)))

#define GET_LAST_5_BITS(ch)        (0x001fU & ch)
#define SET_HIGH5(ch)              (((unsigned short)ch) << 10U)
#define SET_MID5(ch)               (((unsigned short) GET_LAST_5_BITS(ch)) << 5U)
#define PACK_TO_SHORT(ch1, ch2, ch3) (SET_HIGH5(ch1) | SET_MID5(ch2) | GET_LAST_5_BITS(ch3))

#define GET_HIGH_CHAR(num)         (0x40U | ((0x7c00U & num) >> 10U))
#define GET_MID_CHAR(num)          (0x40U | ((0x03e0U & num) >> 5U))
#define GET_LOW_CHAR(num)          (0x40U | (0x001fU & num))

#endif //SRASEARCH_BITMANIPULATEMACROS_H
