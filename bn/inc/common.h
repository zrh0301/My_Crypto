#ifndef __HEADER_COMMON_H__
#define __HEADER_COMMON_H__

#include "wlc_types.h"
/*
*描述：将unsigned char数组转化为dig_t 数
*/
/*
*描述：将32数转化为unsigned char数组
*输入：dig_t数=“87654321 17654321 27654321 37654321”
*输出：unsigned char数组= [0]="87654321"
*               [1]="17654321"
*               [2]="2......."
*               [3]="3......."
*/

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * Convert a character string to hexadecimal data and store it in Byte.
 *
 * @param[in] ch			- characters to be processed.
 * @param[out] ch_byte		- byte type data after processing.
 * @return 1 if the conversion is successful, otherwise 0.
 */
int ConvertHexChar(char ch, unsigned char* ch_byte);

/**
 * Convert a char type array to unsigned char type array.
 *
 * @param[in] pCharBuf		- a char type array.
 * @param[in] charlen		- the length of the char.
 * @param[out] pByteBuf		- a unsigned char type array.
 * @param[out] bytelen		- the length of the unsigned char.
 * @return 1 if the conversion is successful, otherwise 0.
 */
int CharToByte(char* pCharBuf, int charlen, unsigned char* pByteBuf, int* bytelen);

/**
 * Convert unsigned char type array to big integer.
 *
 * @param[in] pByteBuf		- a unsigned char type array.
 * @param[in] bytelen		- the length of the unsigned char.
 * @param[out] pwBN			- a dig_t type array.
 * @param[out] iBNWordLen	- the length of the dig_t.
 * @return 1 if the conversion is successful, otherwise 0.
 */
int ByteToBN(unsigned char* pByteBuf, int bytelen, dig_t* pwBN, int iBNWordLen);

/**
 *change the order of a
 * 
 * @param[out] r			- the result.
 * @param[in] a				- the big number to inverse the order.
 * @param[in] len			- the number of digits.
 */
void inverse(dig_t* r, const dig_t* a, int len);

#ifdef  __cplusplus
}
#endif

#endif
