// WavOverSampling.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>


#include <Windows.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

#if defined(BOOST_VERSION)
namespace mp = boost::multiprecision;
using boost::math::constants::pi;
#endif

// 16bit   scale: 40bit  56bit->32
#define SCALE (1LL << 40)
#define SCALE_SHIFT 24

extern "C" void createHannCoeff(int tapNum, long long* dest, long long scale)
{
	using namespace boost::multiprecision;

	int tapNumMed = ((tapNum - 1) / 2) + 1;
	int coeffNum = (tapNum + 1) / 2;

	scale = SCALE;

	for (int i = 1; i < coeffNum; ++i)
	{

		cpp_dec_float_100 piq = pi<cpp_dec_float_100>();
		cpp_dec_float_100 x = (piq * i) / 4;

		//float128 coeff = (float128)1.0 / (piq * i) * sin(x);
		cpp_dec_float_100 coeff = 1 / (piq * i) * sin(x);

		//float128 hamming = 0.5 + 0.5*cos((i * piq) / (tapNumMed - 1));
		cpp_dec_float_100 a = i;
		cpp_dec_float_100 b = tapNumMed - 1;
		cpp_dec_float_100 hamming("0.5");
		a *= piq;
		a /= b;
		a = cos(a);
		hamming = hamming * a + hamming;

		//float128 coeffAfterWin = coeff * hamming * scale;
		cpp_dec_float_100 coeffAfterWin = scale;
		coeffAfterWin = coeffAfterWin * hamming * coeff;
		coeffAfterWin = round(coeffAfterWin);

		dest[coeffNum - 1 + i] = (long long)coeffAfterWin; // 2 ^ 30
		dest[coeffNum - 1 - i] = (long long)coeffAfterWin;
	}

	dest[coeffNum - 1] = (long long)scale >> 2; ////!!!!

}

static void writeRaw32bitPCM(long long left, long long right, int* buffer)
{
	//left = left >> 12; // 12 is OK for scale 2^30 & 63tap :: <16bit> 30bit - 2bit (max 0.25) - 12 = 32   //// 2 ^ 30
	//0x20000000 30bit 4bit
	// 4095 : 
	int shift = SCALE_SHIFT+1;

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

	left = left >> shift;
	right = right >> shift;


	buffer[0] = (int)left;
	buffer[1] = (int)right;
}


int  oversample(short* src, unsigned int length, long long* coeff, int tapNum, int* dest)
{
	int mid = (tapNum + 1) / 2 - 1; // e.g  15 for 31tap,  63 for 127tap


	for (unsigned int i = 0; i < length; ++i)
	{
		short *srcLeft = src;
		short *srcRight = src + 1;
		long long tmpLeft, tmpRight;

		int startCoeff;
		int n;

		// 1st 
		tmpLeft = *srcLeft * coeff[mid];
		tmpRight = *srcRight * coeff[mid];
		writeRaw32bitPCM(tmpLeft, tmpRight, dest);

		// 2nd  current * [-1], current-1 * [-5], current-2 * [-9]
		tmpLeft = 0;
		tmpRight = 0;
		startCoeff = (mid - 1) % 8;
		n = (tapNum - startCoeff) / 8;
		srcLeft = src - ((mid - 1) / 8) * 2;
		srcRight = srcLeft + 1;

		for (int j = 0; j < n; ++j)
		{
			tmpLeft += coeff[startCoeff + j * 8] * (long long)*(srcLeft + j * 2);
			tmpRight += coeff[startCoeff + j * 8] * (long long)*(srcRight + j * 2);
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 2);

		// 3rd  current * [-1], current-1 * [-5], current-2 * [-9]
		tmpLeft = 0;
		tmpRight = 0;

		startCoeff = (mid - 2) % 8;
		n = (tapNum - startCoeff) / 8;
		srcLeft = src - ((mid - 2) / 8) * 2;
		srcRight = srcLeft + 1;

		for (int j = 0; j < n; ++j)
		{
			tmpLeft += coeff[startCoeff + j * 8] * (long long)*(srcLeft + j * 2);
			tmpRight += coeff[startCoeff + j * 8] * (long long)*(srcRight + j * 2);
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 4);

		// 4th
		tmpLeft = 0;
		tmpRight = 0;
		startCoeff = (mid - 3) % 8;
		n = (tapNum - startCoeff) / 8;
		srcLeft = src - ((mid - 3) / 8) * 2;
		srcRight = srcLeft + 1;
		for (int j = 0; j < n; ++j)
		{
			tmpLeft += coeff[startCoeff + j * 8] * (long long)*(srcLeft + j * 2);
			tmpRight += coeff[startCoeff + j * 8] * (long long)*(srcRight + j * 2);
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 6);

		//5th
		tmpLeft = 0;
		tmpRight = 0;
		startCoeff = (mid - 4) % 8;
		n = (tapNum - startCoeff) / 8;
		srcLeft = src - ((mid - 4) / 8) * 2;
		srcRight = srcLeft + 1;
		for (int j = 0; j < n; ++j)
		{
			tmpLeft += coeff[startCoeff + j * 8] * (long long)*(srcLeft + j * 2);
			tmpRight += coeff[startCoeff + j * 8] * (long long)*(srcRight + j * 2);
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 8);

		//6th
		tmpLeft = 0;
		tmpRight = 0;
		startCoeff = (mid - 5) % 8;
		n = (tapNum - startCoeff) / 8;
		srcLeft = src - ((mid - 5) / 8) * 2;
		srcRight = srcLeft + 1;
		for (int j = 0; j < n; ++j)
		{
			tmpLeft += coeff[startCoeff + j * 8] * (long long)*(srcLeft + j * 2);
			tmpRight += coeff[startCoeff + j * 8] * (long long)*(srcRight + j * 2);
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 10);

		//7th
		tmpLeft = 0;
		tmpRight = 0;
		startCoeff = (mid - 6) % 8;
		n = (tapNum - startCoeff) / 8;
		srcLeft = src - ((mid - 6) / 8) * 2;
		srcRight = srcLeft + 1;
		for (int j = 0; j < n; ++j)
		{
			tmpLeft += coeff[startCoeff + j * 8] * (long long)*(srcLeft + j * 2);
			tmpRight += coeff[startCoeff + j * 8] * (long long)*(srcRight + j * 2);
		}
		writeRaw32bitPCM(tmpLeft, tmpRight, dest + 12);

		//8th
		tmpLeft = 0;
		tmpRight = 0;
		startCoeff = (mid - 7) % 8;
		n = (tapNum - startCoeff) / 8;
		srcLeft = src - ((mid - 7) / 8) * 2;
		srcRight = srcLeft + 1;
		for (int j = 0; j < n; ++j)
		{
			tmpLeft += coeff[startCoeff + j * 8] * (long long)*(srcLeft + j * 2);
			tmpRight += coeff[startCoeff + j * 8] * (long long)*(srcRight + j * 2);
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

void* readWavFile(wchar_t* fileName)
{
	HANDLE fileHandle;
	DWORD offset, size, readSize;
	void* data;
	WAVEFORMATEX wf;

	if (!searchFmtDataChunk(fileName, &wf, &offset, &size))
	{
		return 0;
	}

	fileHandle = CreateFileW(fileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	if (fileHandle == INVALID_HANDLE_VALUE)
	{
		return 0;
	}

	data = (void*) ::GlobalAlloc(GPTR, ((size + 15) / 16) * 16);

	::SetFilePointer(fileHandle, offset, 0, FILE_CURRENT);
	::ReadFile(fileHandle, data, size, &readSize, NULL);
	::CloseHandle(fileHandle);

	return data;
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

int writeWAV_header(HANDLE fileHandle, int ch, int freq, int depth, unsigned long dataSize)
{
	if (fileHandle == 0) return 0;
	if (ch != 1 && ch != 2) return 0;
	if (depth != 16 && depth != 24 && depth != 32) return 0;

	WAVEFORMATEX wf;
	wf.wFormatTag = 0x01; // 1:LPCM
	wf.nChannels = ch;
	wf.nSamplesPerSec = freq;
	wf.nAvgBytesPerSec = freq * ((depth * ch) / 8);
	wf.nBlockAlign = (depth * ch) / 8; // 4bytes (16bit, 2ch) per sample
	wf.wBitsPerSample = depth;
	wf.cbSize = 0;

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


// X次 ノイズシェーパー 
#define DEPTH_SHIFT 8
void noiseShaper_24bitTo16bit(wchar_t* destFileName, void* data, int sample, WAVEFORMATEX* wf, int x)
{
	HANDLE fileHandle = CreateFileW(destFileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS /*CREATE_NEW*/, FILE_ATTRIBUTE_NORMAL, NULL);
	writeWAV_header(fileHandle, 2, wf->nSamplesPerSec, 16, sample * 2 * 2);

	short* data2 = (short*)::GlobalAlloc(GPTR, sample * 2 * 2);

	BYTE* p = (BYTE*)data;
	int left, right, add = 1 << (DEPTH_SHIFT - 1);
	DWORD writtenSize;

	int lastRight = 0;
	int lastLeft = 0;
	int outRight;
	int outLeft;
	long long sigmaLeft[32];
	long long sigmaRight[32];

	short n = -1;

	for (int i = 0; i < 32; ++i)
	{
		sigmaLeft[i] = sigmaRight[i] = 0;
	}

	for (int i = 0; i < sample; ++i)
	{
		unsigned int tmp;
		tmp = (p[2] << 16) | (p[1] << 8) | p[0];
		if (tmp > 0x00800000)
		{
			left = (0x01000000 - tmp) * -1;
		}
		else left = (int)tmp;

		tmp = (p[5] << 16) | (p[4] << 8) | p[3];
		if (tmp > 0x00800000)
		{
			right = (0x01000000 - tmp) * -1;
		}
		else right = (int)tmp;

		sigmaLeft[0] += left - lastLeft;
		sigmaRight[0] += right - lastRight;
		for (int j = 1; j < x; ++j)
		{
			sigmaLeft[j] += sigmaLeft[j - 1] - lastLeft;
			sigmaRight[j] += sigmaRight[j - 1] - lastRight;
		}

		if (sigmaLeft[x - 1] >= 0)
		{
			outLeft = (sigmaLeft[x - 1] + add) >> DEPTH_SHIFT;
			lastLeft = outLeft << 8;
		}
		else
		{
			outLeft = (sigmaLeft[x - 1] * -1 + add) >> DEPTH_SHIFT;
			lastLeft = outLeft << 8;
			outLeft *= -1;
			lastLeft *= -1;
		}

		///// right
		if (sigmaRight[x - 1] >= 0)
		{
			outRight = (sigmaRight[x - 1] + add) >> DEPTH_SHIFT;
			lastRight = outRight << 8;
		}
		else
		{
			outRight = (sigmaRight[x - 1] * -1 + add) >> DEPTH_SHIFT;
			lastRight = outRight << 8;
			outRight *= -1;
			lastRight *= -1;
		}

		if (outLeft > 32767) outLeft = 32767;
		if (outLeft < -32768) outLeft = -32768;

		if (outRight > 32767) outRight = 32767;
		if (outRight < -32768) outRight = -32768;

		data2[i * 2] = (short)outLeft;
		data2[i * 2 + 1] = (short)outRight;

		p += 6;
	}

	::WriteFile(fileHandle, data2, sample * 2 * 2, &writtenSize, NULL);
	::FlushFileBuffers(fileHandle); // なくても大丈夫そう
	::CloseHandle(fileHandle);

	::GlobalFree(data2);
}

#define DATA_UNIT_SIZE (1024 * 1024)
#define TAP_SIZE 4096

int main()
{
	DWORD wavDataOffset, wavDataSize, readSize = 0;
	WAVEFORMATEX wf;
	wchar_t fileName[] = L"C:\\Test\\1k_174_24.wav";
	//wchar_t fileName[] = L"C:\\Test\\a001.WAV";
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
	createHannCoeff(TAP_SIZE, firCoeff, 32);

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

// プログラムの実行: Ctrl + F5 または [デバッグ] > [デバッグなしで開始] メニュー
// プログラムのデバッグ: F5 または [デバッグ] > [デバッグの開始] メニュー

// 作業を開始するためのヒント: 
//    1. ソリューション エクスプローラー ウィンドウを使用してファイルを追加/管理します 
//   2. チーム エクスプローラー ウィンドウを使用してソース管理に接続します
//   3. 出力ウィンドウを使用して、ビルド出力とその他のメッセージを表示します
//   4. エラー一覧ウィンドウを使用してエラーを表示します
//   5. [プロジェクト] > [新しい項目の追加] と移動して新しいコード ファイルを作成するか、[プロジェクト] > [既存の項目の追加] と移動して既存のコード ファイルをプロジェクトに追加します
//   6. 後ほどこのプロジェクトを再び開く場合、[ファイル] > [開く] > [プロジェクト] と移動して .sln ファイルを選択します
