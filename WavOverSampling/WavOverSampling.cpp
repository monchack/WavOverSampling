/**
 * Copyright 2020 Monchack Audio
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <iostream>

#include <Windows.h>

// reomve comment out below to use Boost
//#include <boost/multiprecision/cpp_dec_float.hpp>

// Tap size; change this number if necessary
#define TAP_SIZE 4095


#define DATA_UNIT_SIZE (1024 * 1024)

// 16(15+1)bit  X  scale: 48(47+1)bit =  63(62+1)bit -> 32bit (31bit shift)
#define COEFF_SCALE 47
#define SCALE_SHIFT 31

#if !defined(BOOST_VERSION)

void createHannCoeff(int tapNum, long long* dest)
{
	int coeffNum = (tapNum + 1) / 2;
	double* coeff1 = (double*)::GlobalAlloc(GPTR, sizeof(double) * coeffNum);
	double* coeff2 = (double*)::GlobalAlloc(GPTR, sizeof(double) * coeffNum);
	double* coeff3 = (double*)::GlobalAlloc(GPTR, sizeof(double) * coeffNum);
	double pi = 3.141592653589793;

	coeff1[0] = 2.0f * (22050.0f / 352800.0f);
	for (int i = 1; i < coeffNum; ++i)
	{
		double x = i * 2.0f * pi * (22050.0f / 352800.0f);
		coeff1[i] = sin(x) / (pi * i);
	}

	for (int i = 0; i < coeffNum; ++i)
	{
		double x = 2.0f * pi * i / (double)(tapNum - 1);
		coeff2[i] = 0.5f + 0.5f * cos(x);
	}
	coeff2[coeffNum - 1] = 0;

	long long scale = 1LL << (COEFF_SCALE + 3);

	for (int i = 0; i < coeffNum; ++i)
	{
		coeff3[i] = round(coeff1[i] * coeff2[i] * scale);
	}

	dest[coeffNum - 1] = (long long)coeff3[0];
	for (int i = 1; i < coeffNum; ++i)
	{
		dest[coeffNum - 1 + i] = (long long)coeff3[i];
		dest[coeffNum - 1 - i] = (long long)coeff3[i];
	}
	::GlobalFree(coeff1);
	::GlobalFree(coeff2);
	::GlobalFree(coeff3);
}

#else

using namespace boost::multiprecision;
using boost::math::constants::pi;

void createHannCoeff(int tapNum, long long* dest)
{
	int coeffNum = (tapNum + 1) / 2;
	cpp_dec_float_100* coeff1 = (cpp_dec_float_100*)::GlobalAlloc(GPTR, sizeof(cpp_dec_float_100) * coeffNum);
	cpp_dec_float_100* coeff2 = (cpp_dec_float_100*)::GlobalAlloc(GPTR, sizeof(cpp_dec_float_100) * coeffNum);
	cpp_dec_float_100* coeff3 = (cpp_dec_float_100*)::GlobalAlloc(GPTR, sizeof(cpp_dec_float_100) * coeffNum);

	cpp_dec_float_100 piq = pi<cpp_dec_float_100>();

	coeff1[0] = 2.0f * (22050.0f / 352800.0f);
	for (int i = 1; i < coeffNum; ++i)
	{
		cpp_dec_float_100 x = i * 2 * piq * 22050 / 352800;
		coeff1[i] = boost::multiprecision::sin(x) / (piq * i);
	}

	for (int i = 0; i < coeffNum; ++i)
	{
		cpp_dec_float_100 x = 2.0 * piq * i / (tapNum - 1);
		coeff2[i] = 0.5 + 0.5 * boost::multiprecision::cos(x);
	}
	coeff2[coeffNum - 1] = 0;

	long long scale = 1LL << (COEFF_SCALE + 3);

	for (int i = 0; i < coeffNum; ++i)
	{
		coeff3[i] = boost::multiprecision::round(coeff1[i] * coeff2[i] * scale);
	}

	dest[coeffNum - 1] = (long long)coeff3[0];
	for (int i = 1; i < coeffNum; ++i)
	{
		dest[coeffNum - 1 + i] = (long long)coeff3[i];
		dest[coeffNum - 1 - i] = (long long)coeff3[i];
	}
	::GlobalFree(coeff1);
	::GlobalFree(coeff2);
	::GlobalFree(coeff3);
}

#endif

static void writeRaw32bitPCM(long long left, long long right, int* buffer)
{
	int shift = SCALE_SHIFT;

	if (left >= 4611686018427387904) left = 4611686018427387904 - 1; // over 63bit : limitted to under [1 << 62]   62bit + 1bit
	if (right >= 4611686018427387904) right = 4611686018427387904 - 1;

	if (left < -4611686018427387904) left = -4611686018427387904;
	if (right < -4611686018427387904) right = -4611686018427387904;

	left = left >> shift;
	right = right >> shift;

	buffer[0] = (int)left;
	buffer[1] = (int)right;
}


int  oversample(short* src, unsigned int length, long long* coeff, int tapNum, int* dest, unsigned int option)
{
	int half_size = (tapNum - 1) / 2;
	if (option == 0) option = 0xffff;

	for (unsigned int i = 0; i < length; ++i)
	{
		short *srcLeft = src;
		short *srcRight = src + 1;
		long long tmpLeft, tmpRight;

		if (option & 0x0001)
		{
			// 1st 
			tmpLeft = *srcLeft * coeff[half_size];
			tmpRight = *srcRight * coeff[half_size];
			writeRaw32bitPCM(tmpLeft, tmpRight, dest);
		}

		if (option & 0x0002)
		{
			// 2nd 
			tmpLeft = 0;
			tmpRight = 0;
			// src[1] * coeff[ 7]  +  src[ 2] * coeff[15]  +  src[ 3] * coeff[ 23]  + ...    
			// src[0] * coeff[-1]  +  src[-1] * coeff[-9]  +  src[-2] * coeff[-17]  + ...
			for (int j = 1; (j * 8 - 1) <= half_size; ++j)
			{
				tmpLeft += (long long)*(srcLeft + j * 2) * coeff[half_size + j * 8 - 1];
				tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 1];
			}
			for (int j = 0; (j * 8 + 1) <= half_size; ++j)
			{
				tmpLeft += (long long)*(srcLeft - j * 2) * coeff[half_size - j * 8 - 1];
				tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 1];
			}
			writeRaw32bitPCM(tmpLeft, tmpRight, dest + 2);
		}

		if (option & 0x0004)
		{
			// 3rd 
			tmpLeft = 0;
			tmpRight = 0;
			// src[1] * coeff[ 6]  +  src[ 2] * coeff[ 14]  +  src[ 3] * coeff[ 22]  + ...    
			// src[0] * coeff[-2]  +  src[-1] * coeff[-10]  +  src[-2] * coeff[-18]  + ...
			for (int j = 1; (j * 8 - 2) <= half_size; ++j)
			{
				tmpLeft += (long long)*(srcLeft + j * 2) * coeff[half_size + j * 8 - 2];
				tmpRight += (long long)*(srcRight + j * 2) * coeff[half_size + j * 8 - 2];
			}
			for (int j = 0; (j * 8 + 2) <= half_size; ++j)
			{
				tmpLeft += (long long)*(srcLeft - j * 2) * coeff[half_size - j * 8 - 2];
				tmpRight += (long long)*(srcRight - j * 2) * coeff[half_size - j * 8 - 2];
			}
			writeRaw32bitPCM(tmpLeft, tmpRight, dest + 4);
		}

		if (option & 0x0008)
		{
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
		}

		if (option & 0x0010)
		{
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
		}

		if (option & 0x0020)
		{
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
		}

		if (option & 0x0040)
		{
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
		}
	
		if (option & 0x0080)
		{
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
		}

		src += 2;
		dest += 8 * 2;
	}

	return 0;
}

struct oversample_info
{
	short* src;
	unsigned int length;
	long long* coeff;
	int tapNum;
	int* dest;
	unsigned int option;
};

DWORD WINAPI ThreadFunc(LPVOID arg)
{
	struct oversample_info* info = (struct oversample_info*)arg;
	oversample(info->src, info->length, info->coeff, info->tapNum, info->dest, info->option);
	return 0;
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

int main()
{
	DWORD wavDataOffset, wavDataSize, writtenSize, length, readSize = 0;
	WAVEFORMATEX wf;
	wchar_t fileName[] = L"C:\\Test\\1k_44_16.wav";
	wchar_t destFileName[] = L"C:\\Test\\out2.WAV";

	ULONGLONG elapsedTime = GetTickCount64();

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

	long long* firCoeff = (long long*)::GlobalAlloc(GPTR, sizeof(long long) * TAP_SIZE);
	createHannCoeff(TAP_SIZE, firCoeff);

	for (int i = 0; i <= part; ++i)
	{
		length = readSize;
		::CopyMemory(mem1, mem2, DATA_UNIT_SIZE);
		::CopyMemory(mem2, mem3, DATA_UNIT_SIZE);
		::SecureZeroMemory(mem3, DATA_UNIT_SIZE);
		if (i != part) readSize = readWavFile(fileName, mem3, DATA_UNIT_SIZE * i, DATA_UNIT_SIZE);
		if (i == 0) continue;
	
		struct oversample_info info[8];
		info[0].src = (short* )mem2;
		info[0].length = length / 4;
		info[0].coeff = firCoeff;
		info[0].tapNum = TAP_SIZE;
		info[0].dest = (int* )memOut;
		info[0].option = 0;
		
		// Single thread
		ThreadFunc((LPVOID)&info[0]);

		/*
		// Multi thread (use code below instead of above)
		HANDLE thread[8];
		DWORD threadId[8];
		for (int j = 0; j < 8; ++j)
		{
			info[j] = info[0];
			info[j].option = 1 << j;
			thread[j] = CreateThread(NULL, 0, ThreadFunc, (LPVOID)&info[j], 0, &threadId[j]);
		}
		WaitForMultipleObjects(8, thread, TRUE, INFINITE);
		*/
		
		::WriteFile(fileOut, memOut, length * 8 * 2, &writtenSize, NULL);

	}
	elapsedTime = GetTickCount64() - elapsedTime;
	std::cout << "WavOverSampling: Completed.   " << (elapsedTime/1000) << "." << (elapsedTime % 1000) <<  " msec  \n";

	::FlushFileBuffers(fileOut);
	::CloseHandle(fileOut);

	::GlobalFree(mem1);
	::GlobalFree(memOut);
	::GlobalFree(firCoeff);
}
