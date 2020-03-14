// WavOverSampling.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>

#include <Windows.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

#if defined(BOOST_VERSION)
namespace mp = boost::multiprecision;
using boost::math::constants::pi;
#endif

// 16bit   scale: 46bit  62bit->32bit
#define SCALE (1LL << 46)
#define SCALE_SHIFT 30


void createHannCoeff(int tapNum, long long* dest)
{
	using namespace boost::multiprecision;

	int tapNumMed = ((tapNum - 1) / 2) + 1;
	int coeffNum = (tapNum + 1) / 2;

	long long scale = SCALE;

	for (int i = 1; i < coeffNum; ++i)
	{
		cpp_dec_float_100 piq = pi<cpp_dec_float_100>();
		cpp_dec_float_100 x = (piq * i) / 8;   // 2pi *  (1/16)

		//float128 coeff = (float128)1.0 / (piq * i) * sin(x);
		cpp_dec_float_100 coeff = 1 / (piq * i) * sin(x);

		//float128 hamming = 0.5 + 0.5*cos((i * piq) / (tapNumMed - 1));
		cpp_dec_float_100 a = i;
		cpp_dec_float_100 b = tapNumMed - 1;
		cpp_dec_float_100 hann("0.5");
		a *= piq;
		a /= b;
		a = cos(a);
		hann = hann * a + hann;

		//float128 coeffAfterWin = coeff * hann * scale;
		cpp_dec_float_100 coeffAfterWin = scale;
		coeffAfterWin = coeffAfterWin * hann * coeff;
		coeffAfterWin = round(coeffAfterWin);

		dest[coeffNum - 1 + i] = (long long)coeffAfterWin; // 2 ^ 30
		dest[coeffNum - 1 - i] = (long long)coeffAfterWin;
	}

	dest[coeffNum - 1] = (long long)scale >> 3; ////!!!!

}

static void writeRaw32bitPCM(long long left, long long right, int* buffer)
{
	//left = left >> 12; // 12 is OK for scale 2^30 & 63tap :: <16bit> 30bit - 2bit (max 0.25) - 12 = 32   //// 2 ^ 30
	//0x20000000 30bit 4bit
	// 4095 : 
	int shift = SCALE_SHIFT;

	/*
	if (left > 0)
	{
		left += (1LL << (shift - 1));
	}
	else if (left < 0)
	{
		left -= (1LL << (shift - 1));
	}
	if (right > 0)
	{
		right += (1LL << (shift - 1));
	}
	else if (right < 0)
	{
		right -= (1LL << (shift - 1));
	}
	*/
	

	left = left >> shift;
	right = right >> shift;


	buffer[0] = (int)left;
	buffer[1] = (int)right;
}


int  oversample(short* src, unsigned int length, long long* coeff, int tapNum, int* dest)
{
	for (unsigned int i = 0; i < length; ++i)
	{
		short *srcLeft = src;
		short *srcRight = src + 1;
		long long tmpLeft, tmpRight;

		int half_size = (tapNum - 1) / 2;

		// 1st 
		tmpLeft = *srcLeft * coeff[half_size];
		tmpRight = *srcRight * coeff[half_size];
		writeRaw32bitPCM(tmpLeft, tmpRight, dest);

		// 2nd 
		tmpLeft = 0;
		tmpRight = 0;
		// src[1] * coeff[ 7]  +  src[ 2] * coeff[15]  +  src[ 3] * coeff[ 23]  + ...    
		// src[0] * coeff[-1]  +  src[-1] * coeff[-9]  +  src[-2] * coeff[-17]  + ...
		for (int j = 1; (j * 8 - 1) <= half_size; ++j)
		{
			tmpLeft  += (long long)*(srcLeft  + j * 2) * coeff[half_size + j * 8 - 1];
			tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 1];
		}
		for (int j = 0; (j * 8 + 1 ) <= half_size; ++j)
		{
			tmpLeft  += (long long)*(srcLeft  - j * 2) * coeff[half_size - j * 8 - 1];
			tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 1];
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 2);

		// 3rd 
		tmpLeft  = 0;
		tmpRight = 0;
		// src[1] * coeff[ 6]  +  src[ 2] * coeff[ 14]  +  src[ 3] * coeff[ 22]  + ...    
		// src[0] * coeff[-2]  +  src[-1] * coeff[-10]  +  src[-2] * coeff[-18]  + ...
		for (int j = 1; (j * 8 - 2) <= half_size; ++j)
		{
			tmpLeft  += (long long)*(srcLeft  + j * 2) * coeff[half_size + j * 8 - 2];
			tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 2];
		}
		for (int j = 0; (j * 8 + 2) <= half_size; ++j)
		{
			tmpLeft  += (long long)*(srcLeft  - j * 2) * coeff[half_size - j * 8 - 2];
			tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 2];
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 4);

		// 4th
		tmpLeft = 0;
		tmpRight = 0;
		for (int j = 1; (j * 8 - 3) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft + j * 2) * coeff[half_size + j * 8 - 3];
			tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 3];
		}
		for (int j = 0; (j * 8 + 3) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft - j * 2) * coeff[half_size - j * 8 - 3];
			tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 3];
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 6);

		//5th
		tmpLeft = 0;
		tmpRight = 0;
		for (int j = 1; (j * 8 - 4) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft + j * 2) * coeff[half_size + j * 8 - 4];
			tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 4];
		}
		for (int j = 0; (j * 8 + 4) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft - j * 2) * coeff[half_size - j * 8 - 4];
			tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 4];
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 8);

		//6th
		tmpLeft = 0;
		tmpRight = 0;
		for (int j = 1; (j * 8 - 5) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft + j * 2) * coeff[half_size + j * 8 - 5];
			tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 5];
		}
		for (int j = 0; (j * 8 + 5) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft - j * 2) * coeff[half_size - j * 8 - 5];
			tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 5];
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 10);

		//7th
		tmpLeft = 0;
		tmpRight = 0;
		for (int j = 1; (j * 8 - 6) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft + j * 2) * coeff[half_size + j * 8 - 6];
			tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 6];
		}
		for (int j = 0; (j * 8 + 6) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft - j * 2) * coeff[half_size - j * 8 - 6];
			tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 6];
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 12);

		//8th
		tmpLeft = 0;
		tmpRight = 0;
		for (int j = 1; (j * 8 - 7) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft + j * 2) * coeff[half_size + j * 8 - 7];
			tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 7];
		}
		for (int j = 0; (j * 8 + 7) <= half_size; ++j)
		{
			tmpLeft += (long long)*(srcLeft - j * 2) * coeff[half_size - j * 8 - 7];
			tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 7];
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 14);

		src += 2;
		dest += 8 * 2;
	}

	return 0;
}

int getFileSize(wchar_t* fileName, DWORD* sizeLow, DWORD* sizeHigh)
{
	HANDLE h;
	DWORD dwSizeHigh;
	DWORD dwSizeLow;
	DWORD dwError;

	if (sizeLow == 0 || sizeHigh == 0) return 0;
	*sizeLow = *sizeHigh = 0;
	h = CreateFileW(fileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	if (h == INVALID_HANDLE_VALUE) return 0;
	dwSizeLow = GetFileSize(h, &dwSizeHigh);
	dwError = GetLastError();
	CloseHandle(h);
	if (dwSizeLow == 0xffffffff && dwError != NO_ERROR) return 0;
	*sizeLow = dwSizeLow;
	*sizeHigh = dwSizeHigh;
	return 1;
}

unsigned int searchFmtDataChunk(wchar_t* fileName, WAVEFORMATEX* wf, DWORD* offset, DWORD* size)
{
	HANDLE fileHandle;
	fileHandle = CreateFileW(fileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	if (fileHandle == INVALID_HANDLE_VALUE)
	{
		return 0;
	}

	DWORD header[2];
	DWORD readSize;
	WORD  wav[8];
	DWORD riffSize, pos = 0;
	DWORD dataOffset, dataSize;
	::ReadFile(fileHandle, header, 8, &readSize, NULL);
	bool fmtFound = false, dataFound = false;

	if (readSize != 8)
	{
		CloseHandle(fileHandle);
		return 0;
	}

	if (header[0] != 0X46464952)
	{
		// not "RIFF"
		CloseHandle(fileHandle);
		return 0;
	}
	riffSize = header[1];

	::ReadFile(fileHandle, header, 4, &readSize, NULL);
	if (readSize != 4)
	{
		CloseHandle(fileHandle);
		return 0;
	}
	if (header[0] != 0x45564157)
	{
		// not "WAVE"
		CloseHandle(fileHandle);
		return 0;
	}
	pos += 4;

	while (pos < riffSize)
	{
		::ReadFile(fileHandle, header, 8, &readSize, NULL);
		if (readSize != 8)
		{
			break;
		}
		pos += 8;

		if (header[0] == 0X20746d66)
		{
			// "fmt "
			if (header[1] >= 16)
			{
				::ReadFile(fileHandle, wav, 16, &readSize, NULL);
				if (readSize != 16)
				{
					break;
				}
				fmtFound = true;
				if (header[1] > 16)
				{
					::SetFilePointer(fileHandle, header[1] - 16, 0, FILE_CURRENT);
				}
				pos += header[1];
			}
			else
			{
				::SetFilePointer(fileHandle, header[1], 0, FILE_CURRENT);
				pos += header[1];
			}
		}
		else if (header[0] == 0X61746164)
		{
			// "data"
			dataFound = true;
			dataOffset = ::SetFilePointer(fileHandle, 0, 0, FILE_CURRENT);
			dataSize = header[1];
			::SetFilePointer(fileHandle, header[1], 0, FILE_CURRENT);
			pos += header[1];
		}
		else
		{
			::SetFilePointer(fileHandle, header[1], 0, FILE_CURRENT);
			pos += header[1];
		}
		if (GetLastError() != NO_ERROR)
		{
			break;
		}
	}
	CloseHandle(fileHandle);

	if (dataFound && fmtFound)
	{
		*offset = dataOffset;
		*size = dataSize;
		wf->wFormatTag = wav[0]; //  1:LPCM   3:IEEE float
		wf->nChannels = wav[1]; //  1:Mono  2:Stereo
		wf->nSamplesPerSec = *(DWORD*)(wav + 2);  // 44100, 48000, 176400, 19200, 352800, 384000...
		wf->nAvgBytesPerSec = *(DWORD*)(wav + 4);
		wf->nBlockAlign = wav[6]; // 4@16bit/2ch,  6@24bit/2ch,   8@32bit/2ch   
		wf->wBitsPerSample = wav[7]; // 16bit, 24bit, 32bit
		wf->cbSize = 0;
		return 1;
	}
	return 0;
}

DWORD readWavFile(wchar_t* fileName, void* readMem, DWORD readPos, DWORD readLength)
{
	HANDLE fileHandle;
	DWORD wavDataOffset, wavDataSize, readSize = 0;
	WAVEFORMATEX wf;

	if (!searchFmtDataChunk(fileName, &wf, &wavDataOffset, &wavDataSize))
	{
		return 0;
	}

	fileHandle = CreateFileW(fileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	if (fileHandle == INVALID_HANDLE_VALUE)
	{
		return 0;
	}

	if (::SetFilePointer(fileHandle, wavDataOffset + readPos, 0, FILE_BEGIN) == INVALID_SET_FILE_POINTER)
	{
		if (GetLastError() != NO_ERROR)
		{
			// fail
			return 0;
		}
	}
	::ReadFile(fileHandle, readMem, readLength, &readSize, NULL);
	::CloseHandle(fileHandle);

	return readSize;
}

static int writePCM352_32_header(HANDLE fileHandle, unsigned long dataSize)
{
	WAVEFORMATEX wf;
	wf.wFormatTag = 0x01;
	wf.nChannels = 2;
	wf.nSamplesPerSec = 352800;
	wf.nAvgBytesPerSec = 352800 * 8; // 352800 * 4byte(32bit) * 2ch
	wf.nBlockAlign = 8; // 8bytes (32bit, 2ch) per sample
	wf.wBitsPerSample = 32;
	wf.cbSize = 0; // ignored. not written.

	DWORD writtenSize = 0;
	WriteFile(fileHandle, "RIFF", 4, &writtenSize, NULL);
	DWORD size = (dataSize + 44) - 8;
	WriteFile(fileHandle, &size, 4, &writtenSize, NULL);
	WriteFile(fileHandle, "WAVE", 4, &writtenSize, NULL);
	WriteFile(fileHandle, "fmt ", 4, &writtenSize, NULL);
	size = 16;
	WriteFile(fileHandle, &size, 4, &writtenSize, NULL);
	WriteFile(fileHandle, &wf, size, &writtenSize, NULL);
	WriteFile(fileHandle, "data", 4, &writtenSize, NULL);
	size = (DWORD)dataSize;
	WriteFile(fileHandle, &size, 4, &writtenSize, NULL);

	return 0;
}



#define DATA_UNIT_SIZE (1024 * 1024)
#define TAP_SIZE 4095

int main()
{
	DWORD wavDataOffset, wavDataSize, readSize = 0;
	WAVEFORMATEX wf;
	wchar_t fileName[] = L"C:\\Test\\1k_44_16.wav";
	//wchar_t fileName[] = L"C:\\Test\\noise44.WAV";
	wchar_t destFileName[] = L"C:\\Test\\out2.WAV";
	DWORD writtenSize;

	if (!searchFmtDataChunk(fileName, &wf, &wavDataOffset, &wavDataSize))
	{
		return 0;
	}
	int part = wavDataSize / DATA_UNIT_SIZE;
	if ((wavDataSize %  DATA_UNIT_SIZE) != 0) part += 1;

	void* mem1 = ::GlobalAlloc(GPTR, DATA_UNIT_SIZE * 3);
	void* mem2 = (char*)mem1 + DATA_UNIT_SIZE;
	void* mem3 = (char*)mem2 + DATA_UNIT_SIZE;

	void* memOut = ::GlobalAlloc(GPTR, DATA_UNIT_SIZE * 8 * 2);

	HANDLE fileOut = CreateFileW(destFileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS /*CREATE_NEW*/, FILE_ATTRIBUTE_NORMAL, NULL);
	writePCM352_32_header(fileOut, wavDataSize * 8 * 2);

	long long firCoeff[TAP_SIZE];
	createHannCoeff(TAP_SIZE, firCoeff);

	for (int i = 0; i < part; ++i)
	{
		DWORD readSize;
		if (i != 0)        readWavFile(fileName, mem1, DATA_UNIT_SIZE * (i - 1), DATA_UNIT_SIZE);
		readSize =         readWavFile(fileName, mem2, DATA_UNIT_SIZE * i,       DATA_UNIT_SIZE);
		if (i != part - 1) readWavFile(fileName, mem3, DATA_UNIT_SIZE * (i + 1), DATA_UNIT_SIZE);

		//int  oversample(short* src, unsigned int length, long long* coeff, int tapNum, int* dest)
		oversample((short* )mem2, DATA_UNIT_SIZE / 4, firCoeff, TAP_SIZE, (int* )memOut);
		::WriteFile(fileOut, memOut, readSize * 8 * 2, &writtenSize, NULL);

	}
	std::cout << "Hello World!\n";

	::FlushFileBuffers(fileOut); // なくても大丈夫そう
	::CloseHandle(fileOut);

	//	HANDLE fileHandle = CreateFileW(destFileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS /*CREATE_NEW*/, FILE_ATTRIBUTE_NORMAL, NULL);

	::GlobalFree(mem1);
}
