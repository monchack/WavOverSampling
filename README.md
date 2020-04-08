# WavOverSampling
Oversample 44.1kHz WAV file to 352.8kHz

Tap number is variable but givein in the source code. It must be an odd number.
Also high precision mode can be enabled easily by modifying the code.
Input wav file must be 44.1/16bit/Stereo audio file.
Output wav file's audio format is fixed; 352.8/32bit/Stereo.

usage:
WavOverSampling.exe  inputfile.wav  outputfile.wav
