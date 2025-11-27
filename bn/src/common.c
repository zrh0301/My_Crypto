#include "../inc/common.h"

#define LSBOfWord					0x00000001
#define MSBOfWord					0x80000000

/**
 * Convert a character string to hexadecimal data and store it in Byte.
 *
 * @param[in] ch			- characters to be processed.
 * @param[out] ch_byte		- byte type data after processing.
 * @return 1 if the conversion is successful, otherwise 0.
 */
int ConvertHexChar(char ch, unsigned char* ch_byte)
{
	if ((ch >= '0') && (ch <= '9'))
	{

		*ch_byte = (unsigned char)(ch - 0x30);
		return 1;

	}
	else
	{
		if ((ch >= 'A') && (ch <= 'F'))
		{
			*ch_byte = (unsigned char)(ch - 'A' + 0x0a);
			return 1;
		}
		else
		{
			if ((ch >= 'a') && (ch <= 'f'))
			{
				*ch_byte = (unsigned char)(ch - 'a' + 0x0a);
				return 1;
			}
		}
	}
	return 0;
}

/**
 * Convert a char type array to unsigned char type array.
 *
 * @param[in] pCharBuf		- a char type array.
 * @param[in] charlen		- the length of the char.
 * @param[out] pByteBuf		- a char type array.
 * @param[out] bytelen		- the length of the unsigned char.
 * @return 1 if the conversion is successful, otherwise 0.
 */
int CharToByte(char *pCharBuf, int charlen, unsigned char *pByteBuf, int *bytelen)
{
	/*******************/
	int i = 0;
	int charlen_tmp = 0;
	unsigned char hdata = 0;
	unsigned char ldata = 0;
	char hstr = 0;
	char lstr = 0;
	/*******************/


	charlen_tmp = charlen;
	if (charlen_tmp & LSBOfWord)
	{
		charlen_tmp += 1;
		*bytelen = charlen_tmp >> 1;
		if (ConvertHexChar(pCharBuf[0], &ldata) == 1)
		{
			pByteBuf[0] = ldata;
		}
		for (i = 1; i < *bytelen; i++)
		{
			if (ConvertHexChar(pCharBuf[2 * i - 1], &hdata) == 1)
			{
				if (ConvertHexChar(pCharBuf[2 * i], &ldata) == 1)
				{
					pByteBuf[i] = (hdata << 4) | ldata;
				}
				else
				{
					return 0;
				}
			}
			else
			{
				return 0;
			}
		}		 
	}
	else
	{
		*bytelen = charlen_tmp >> 1;
		for (i = 0; i < *bytelen; i++)
		{
			if (ConvertHexChar(pCharBuf[2 * i], &hdata) == 1)
			{
				if (ConvertHexChar(pCharBuf[2 * i + 1], &ldata) == 1)
				{
					pByteBuf[i] = (hdata << 4) | ldata;
				}
				else
				{
					return 0;
				}
			}
			else
			{
				return 0;
			}
		}
	}	
	return 1;
}

/**
 * Convert unsigned char type array to big integer.
 *
 * @param[in] pByteBuf		- a unsigned char type array.
 * @param[in] bytelen		- the length of the unsigned char.
 * @param[out] pwBN			- a dig_t type array.
 * @param[out] iBNWordLen	- the length of the dig_t.
 * @return 1 if the conversion is successful, otherwise 0.
 */
int ByteToBN(unsigned char *pByteBuf, int bytelen, dig_t *pwBN, int iBNWordLen)
{
	/*******************/
	int ExpLen = 0;
	int Rem = 0;
	int i = 0;
	int j = 0; 
	/*******************/

	ExpLen = bytelen >> 2;
	Rem = bytelen & 0x00000003;

	if (Rem != 0)
	{
		ExpLen += 1; 
	}

	if (ExpLen > iBNWordLen)
	{
		return 0;
	}

	i = bytelen - 1;
	j = 0;
	while (i >= Rem)
	{
		pwBN[j] = ((dig_t)pByteBuf[i]) | ((dig_t)pByteBuf[i - 1] << 8) | ((dig_t)pByteBuf[i - 2] << 16) | ((dig_t)pByteBuf[i - 3] << 24);
		i -= 4;
		j++;
	}

	i = 0;
	while (i < Rem)
	{
		pwBN[j] = (pwBN[j] << 8) | ((dig_t)pByteBuf[i]);
		i++;
	}

	return 1;
}

/**
 *change the order of a
 *
 * @param[out] r			- the result.
 * @param[in] a				- the big number to inverse the order.
 * @param[in] len			- the number of digits.
 */
void inverse(dig_t* r, const dig_t* a, int len)
{
	int temp;
	for (int i = 0; i < len / 2; i++)
	{
		temp = a[i];
		r[i] = a[len - i - 1];
		r[len - i - 1] = temp;
	}
}
